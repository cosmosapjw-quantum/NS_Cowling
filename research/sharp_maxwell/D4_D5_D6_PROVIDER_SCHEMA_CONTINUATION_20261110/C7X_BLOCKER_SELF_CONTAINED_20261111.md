# C7X Sharp-Maxwell source-admission blocker: standalone technical record

**Snapshot:** evidence gathered 2026-11-10, consolidated 2026-11-11  
**Disposition:** STOP_BEFORE_D6_ROOT_SEARCH  
**Scope:** CompOSE source admission for a cold Sharp-Maxwell crust-core transition, before phase-aware TOV/Cowling calculations.

## 1. Executive decision

The current blocker is upstream of the sharp-interface stellar solver. The selected source tables can be parsed, but their composition fields cannot yet be assigned a provider-confirmed physical meaning consistent with the standard CompOSE schema and the native verifier's closure equations. Independent thermodynamic and interpolation gates also remain failed or incomplete. The 55-pair filtered D6 screen therefore stops before pressure-equality evaluation and before every Maxwell root search.

This is a data semantics, provenance, and numerical-admission blocker. It is **not** evidence that the microscopic EOSs are physically unstable, and it is **not** evidence that the phase-aware TOV/Cowling solver is invalid. No raw EOS table or production mapping is repaired in this record.

## 2. Physical target and admission condition

For a cold one-conserved-charge Maxwell transition between low-density (−) and high-density (+) phases, coexistence requires

\[
P_-=P_+=P_t,\qquad \mu_{B,-}=\mu_{B,+}=\mu_t.
\]

The jump is one-sided and must satisfy

\[
n_{B,+}>n_{B,-},\qquad
\Delta\epsilon=\epsilon_+-\epsilon_-=\mu_t(n_{B,+}-n_{B,-})=\mu_t\Delta n_B.
\]

A nonempty intersection of tabulated chemical-potential supports is only a precondition. It is not a common-pressure root. Each candidate pair must be evaluated at the same physical chemical potential after documented rest-mass rebasing, then tested for pressure equality, a positive density jump, and the first-law jump identity. Only an admitted pair can be sent to a phase-labelled TOV/Cowling calculation.

## 3. CompOSE group-composition contract

For each particle species CompOSE defines the number fraction \(Y_i=n_i/n_b\), with baryon normalization \(\sum_i B_iY_i=1\). Its mass fraction is a distinct quantity, \(X_i=B_iY_i\).

The nuclear-group record is \((I,A_{\rm av},Z_{\rm av},Y_I)\), where

\[
Y_I=\sum_{j\in M_I}Y_j
=\frac{\sum_j A_jY_j}{A_{\rm av}}.
\]

Thus the group contributes \(A_{\rm av}Y_I\) to baryon normalization and \(Z_{\rm av}Y_I\) to charge. Its mass fraction is \(X_I=A_{\rm av}Y_I\), not \(Y_I\). The provider may define model-specific group indices and phase labels; no retrieved provider sheet declares a different normalization of the fourth group field.

For the inspected four particle pairs, IDs are 10=n, 11=p, 0=e, 1=mu. The standard row closures are

\[
R_B=Y_n+Y_p+A_{\rm av}Y_I-1,\qquad
R_Q=Y_p+Z_{\rm av}Y_I-Y_e-Y_\mu.
\]

## 4. Directly observed D4 evidence

Thirteen eos.compo files were inspected: eleven GPPVA/SPG crust sources plus RG_SLY4 and RG_KDE0v1.

- 4,975 rows total; 4,975 parseable; zero malformed rows.
- Zero inb sequence mismatches against aligned eos.thermo keys.
- The eleven crust tables contain 3,538 rows: 2,387 outer-crust Iphase=0 rows and 1,151 inner-crust Iphase=7 rows.
- In all 1,151 inner-crust rows, the final nuclear-group field equals Y_e exactly at printed precision; Y_n=Y_p=0.
- Neither the standard CompOSE number-fraction closure nor the simple mass-fraction reinterpretation closes both baryon and charge constraints on those 1,151 rows.
- Outer-crust GPPVA/SPG rows instead have final field equal to 1, with Y_e approximately Z_av/A_av. A mass-fraction-like diagnostic closes to printed precision, but this does not establish an alternate provider convention.
- RG(SLy4,KDE0v1) show a different, mostly mass-fraction-like diagnostic pattern. RG_SLy4 heavy-nucleus terminal row 143 is an outlier with r_B=-0.06445 and r_Q=-0.0226704969. RG_KDE0v1 q_extra is 0 for 1,220 rows and 1 for 20 rows, with a transition between rows 1220 and 1221; its meaning is undocumented.

The observations indicate a schema/export mismatch or missing provider documentation, but do not identify the software stage responsible. No phase-dependent conversion is authorized.

## 5. Native verifier evidence and the C5/C6/C6B distinction

The native verifier treats its internal Y_h as a nuclear-group baryon/mass fraction:

\[
r_{C3}=Y_p+(Z_{\rm av}/A_{\rm av})Y_h-Y_e-Y_\mu,\qquad
r_{C4}=Y_n+Y_p+Y_h-1.
\]

C3 is charge closure; C4 is baryon normalization. The C5 implementation mapped Y_h=Y_I for serialization, but exact reader-writer-reader equality proved only round-trip preservation, not physical semantics. The earlier C5 assertion that CompOSE Eq. (4.22) made Y_I itself a mass fraction was incorrect and is superseded. Blindly applying A_av Y_I is not a global repair either: for outer-crust rows with (A_av,Y_I)=(56,1), it gives 56.

A distinct C6 execution used gfortran 13.3.0 and Python 3.12.3. Native build passed, the committed test file passed 6/6, and reader/writer/re-reader/verifier execution completed for SPG_M1 and GPPVA_FSU2. C3 and C4 failed on every Iphase=7 row: 148/148 for SPG_M1 and 52/52 for GPPVA_FSU2. Every Iphase=0 row passed. Source and legacy-output hashes were unchanged.

A subsequent read-only C6B adjudication pinned and independently parsed eos.nb, eos.thermo, eos.compo, eos.t and eos.yq for those two sources. It verified the CompOSE Eq. (4.19) 18-field layout and row keys (iT,inb,iYq) against eos.thermo, rejecting parser/index-alignment as the explanation for those two files. This two-source review is deeper evidence for SPG_M1 and GPPVA_FSU2 only; it is not a replacement for the broader 13-file D4 audit and must not be generalized to other sources without evidence.

C6B details:
- In rows 1-217 of both files, the final group slot is 1, Y_n=Y_p=0 and Y_e=Z/A. Direct Y_h=Y_I closes the native checks there, while the standard group number fraction for a single group would be 1/A and its mass fraction A Y_I.
- Across 200/200 inner-crust rows in those two sources, the stored final slot equals Y_e+Y_mu at printed precision; Y_n=Y_p=0 in all 634 rows examined.
- GPPVA_FSU2 inner-crust rows have integer A,Z and Y_e=Z/A within 5e-6. C3- and C4-implied Y_h agree at about 0.9996-1.0004 in 52/52 rows, while the stored slot is Y_e rather than that value.
- SPG_M1 inner-crust rows have noninteger A,Z and Z/A differs from Y_e. C3 implies Y_h around 0.096-0.543, whereas C4 requires Y_h=1 in all 148 rows. The recorded fields leave about 0.46-0.90 of baryons unaccounted for.
- SPG(M1)'s global parameters match Scurto, Pais and Gulminelli, PRD 109, 103015 (2024), arXiv:2402.15548, whose compressible-liquid-drop cell has cluster and gas components. This is a provenance clue, not a documented definition of the exported group field. The inspected GPPVA literature did not provide rowwise composition values for comparison.

C6B ended PROVIDER_SEMANTICS_UNRESOLVED__BLOCKED. No Y_h=1 rule, reconstructed gas fraction, or provider-specific phase switch is approved.

## 6. Provider documentation boundary

The inspected GPPVA FSU2H sheet identifies Iphase=0/7/6 as outer/inner/homogeneous and separately documents an eos.thermo Nadd geometry integer: 0 outer, 1 droplets, 2 rods, 3 slabs, 4 tubes, 5 bubbles, 6 homogeneous core. Nadd is not the fourth eos.compo field. The SPG M2 sheet documents four particle pairs, one heavy-nucleus quadruple and phase 0/7. The RG_SLy4 sheet documents phase 4 as heavy-nucleus-present and phase 3 as homogeneous. None of the inspected sheets redefines the fourth nuclear-group field. FSU2H's Nadd mapping is not automatically established for every GPPVA sibling, and RG_KDE0v1 q_extra remains opaque.

## 7. D3 thermodynamic admission

The raw-table reconstruction uses

\[
P=n_BQ_1,\quad
\epsilon=n_Bm_nc^2(1+Q_7),\quad
\mu_{Q3}=m_nc^2(1+Q_3),
\]
\[
\mu_E=(\epsilon+P)/n_B=m_nc^2(1+Q_7)+Q_1.
\]

The strict consistency gate is

\[
\frac{|\mu_{Q3}-\mu_E|}
{\max(|\mu_{Q3}|,|\mu_E|,m_nc^2)}\le 2\times10^{-6}.
\]

The CompOSE columns are Q1=P/n_B, Q2=s/n_B, Q3=mu_B/(m_n c^2)-1, Q4=mu_q/(m_n c^2), Q5=mu_l/(m_n c^2), Q6=f/(n_B m_n c^2)-1, and Q7=epsilon/(n_B m_n c^2)-1. The cold neutrino-transparent beta-equilibrium path requires T=0, Q2 approximately 0, Q5 approximately 0 and Q6 approximately Q7. Nonzero Q5 may represent a different lepton/neutrino state and is not silently projected to zero.

In the 16-source C7X subset, strict Q3/Euler fails for seven sources: GPPVA_DD2, GPPVA_DDME2, GPPVA_FSU2R, GPPVA_NL3wrL55, GPPVA_TM1e, RG_SLY4 and RG_KDE0v1. Q5 fails for RG_SLY4 and RG_KDE0v1. A distinct density-reconstruction drift of approximately 7.3e-4 to 7.7e-4 was reported by earlier direct-provider receipts; its admission tolerance remains unresolved pending refinement and observable-sensitivity tests.

## 8. D5 topology-aware Hermite audit

The earlier Hermite replay split only on Iphase. Seven GPPVA tables have additional discrete Nadd geometry boundaries; RG_KDE0v1 has one q_extra 0-to-1 boundary. Including these with phase boundaries gives 45 segments and 15 excluded boundary cells. Seven of those 15 boundary cells would fail the selected interpolation test and eight would pass.

Updated aggregate: 5,535 topology-preserving cells; 1,153 causal-polynomial violations; 1,110 cells with dn/dmu_E <= 0; zero nonpositive interpolated-density cells; 4,382 pass the three selected interval conditions.

**Qualification:** this is an exact boundary-delta correction. The 15 newly discovered boundary cells were individually evaluated; the remaining classifications were inherited from the prior 5,550-cell baseline. It is not a fresh full replay of all retained source cells on a local raw-source tree. The updated script passed syntax/synthetic checks but was not run against a local raw-source project checkout in that session. All eleven crust inputs fail the selected topology-aware D5 exploratory interval gate. These interpolation-polynomial failures are not claims of microscopic EOS instability.

## 9. D6 filtered pair screen and stop rule

The source manifest contains 65 candidates (13 crusts × 5 cores). The current filtered D6 precondition screen contains 55 candidate pairs. In that 55-pair subset:

- 55/55 have nonempty tabulated chemical-potential support overlap.
- 55/55 require rest-mass rebasing.
- 18/55 have both sides passing Q3 and Q5.
- 33/55 have both sides passing Q5.
- 0/55 have both sides passing topology-aware D5.
- All 55 retain a crust-side D4 composition/provider blocker.

Support overlap is not a common-pressure root bracket. No common-mu pressure evaluation, density drift, transition pressure, source-backed density/energy jump, or Maxwell root search was performed. The final decision is STOP_BEFORE_D6_ROOT_SEARCH. D7 segment-preserving adapter, D8 one-pair pilot, D9 cohort sweep and D10 source-backed TOV/Cowling were not run.

## 10. Why no speculative repair is allowed

The current records admit several hypotheses: provider-specific abundance normalization, Wigner-Seitz cell normalization, missing cluster/gas components, provider-side export defect, or a field association not captured by the current data contract. The evidence does not distinguish these uniquely. A remapping chosen only because it reduces closure residuals would silently alter physical composition and potentially change inferred density, pressure, coexistence location and downstream stellar observables.

Keep the raw composition quarantined from production adapters. Do not edit raw EOS, relax verifier tolerances, impose a phase-dependent formula, or initiate roots until the missing semantics and upstream gates are resolved.

## 11. Exact resolution gates

D4 requires provider/CompOSE documentation of the exact fourth-field meaning, denominator and units by model/phase, plus baryon and charge closure equations; alternatively, corrected versioned source tables with provenance and reproducible standard closure. Clarification must address the free n/p denominators, the repeated Y_I=Y_e inner-crust pattern, the SPG_M1 unaccounted gas/baryon component, RG_SLy4 terminal row 143, all seven GPPVA Nadd mappings and RG_KDE0v1 q_extra.

Before root search, also require: (i) D3 Q3/Euler and Q5 failures resolved or explicitly quarantined with justification; (ii) full one-shot topology-preserving D5 replay over every retained source cell with exact source/script hashes; (iii) native parser/read/build on the exact inputs with toolchain provenance; and (iv) D6 rest-mass rebasing, common-mu pressure comparison and density-jump preconditions. Only after a validated pair root may the phase-chain adapter, one-pair pilot, cohort sweep, stable TOV sequence and separate slow/rapid Cowling spectra proceed.

No provider clarification request has been sent in the recorded work.

## 12. Provenance, evidence states and references

- Repository: cosmosapjw-quantum/NS_Cowling.
- Source-authority branch: research/maxwell-native-fortran-v1, commit a09ec30d74b363652d2d1aefe990557c468679dd.
- C7X parent evidence branch: research/sharp-maxwell-independent-source-audit-v1, commit 3dd524e317fde40a07a954e70bc04a4157807f5b.
- Evidence-only continuation branch: research/sharp-maxwell-d4-provider-semantics-v1; PR #5 is open, draft and unmerged.
- C6 native execution identity: 43eea8bd9bee7ebb222e3c66babf5078b398d7ae.
- C6B two-source audit identity: de029388e78f03f9da2aceef78621e7947e1a299.
- Raw EOS tables are not duplicated in this document or its helper bundle; their source-authority identities and row receipts belong to the existing evidence package.

Evidence classes: raw counts and printed-value equality are directly observed; CompOSE meanings are literature-supported by Typel et al., EPJ A 58, 221 (2022), DOI 10.1140/epja/s10050-022-00847-y; closure residuals are derived from the printed rows; C6 build/test/round-trip are execution-verified; D5 is delta-corrected exploratory; D6 roots and D7-D10 are not run.

Primary references:
- CompOSE Reference Manual v3.01, DOI 10.1140/epja/s10050-022-00847-y: https://doi.org/10.1140/epja/s10050-022-00847-y
- CompOSE provider guide: https://compose.obspm.fr/download/pdf/CompOSE_Quick_Guide_for_Providers.pdf
- GPPVA FSU2H data sheet: https://compose.obspm.fr/download/1D/Crust/FSU2H/eos.pdf
- SPG M2 data sheet: https://compose.obspm.fr/download/1D/Crust/M2/eos.pdf
- SPG M1 data sheet: https://compose.obspm.fr/download/1D/Crust/M1/eos.pdf
- RG(SLy4) model page: https://compose.obspm.fr/eos/134
- Scurto, Pais and Gulminelli (2024): https://arxiv.org/abs/2402.15548
- C7X D4-D6 evidence PR #5: https://github.com/cosmosapjw-quantum/NS_Cowling/pull/5

**Final disposition:** preserve the raw data and contradictions; do not infer an undocumented composition mapping; do not run D6 roots.