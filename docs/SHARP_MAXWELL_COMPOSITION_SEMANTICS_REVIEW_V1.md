# Sharp Maxwell Composition Semantics Review (v1)

Status: **OPEN / BLOCKED_REVIEW_REQUIRED**  
Selected sources: \`SPG_M1\`, \`GPPVA_FSU2\` only  
Purpose: adjudicate the meaning and validity of the \`eos.compo\` nuclear-group abundance fields before changing the native composition mapping or entering C7.

## 1. Decision summary

C6 demonstrated that the native Fortran build, committed regression tests, read/write/read path and verifier invocation work for both selected sources. It also demonstrated that C3 and C4 fail on every \(Iphase=7\) row and pass on every \(Iphase=0\) row.

The C5 claim that CompOSE Eq. (4.22) defines \(Y_I\) itself as a baryon/mass fraction was incorrect. Under the standard CompOSE convention, \(Y_I\) is the combined **number fraction** of a nuclear group, and the corresponding baryon/mass fraction is \(X_h=A_{\rm av}Y_I\). However, blindly applying that conversion to the selected raw files yields values greater than one even in outer-crust rows. The current C5 mapping \(Y_h=Y_I\) reproduces the raw field in a round-trip but fails C3/C4 on all inner-crust rows.

Therefore, neither a global direct mapping nor a global \(A_{\rm av}\) conversion is presently justified. The remaining question concerns provider output semantics, possible source-table inconsistency or key/field interpretation, and/or the applicability of the verifier's composition formula to these provider records. No code or physical definition is to be changed until that question is resolved from provider-level evidence.

## 2. Authority and definitions

### 2.1 Standard CompOSE definitions

Primary authority: [CompOSE Reference Manual v3.01](https://doi.org/10.1140/epja/s10050-022-00847-y).

- Eq. (3.10): \(Y_i=n_i/n_b\), particle number fraction.
- Eq. (3.11): \(\sum_i B_iY_i=1\), baryon-number normalization.
- Eq. (3.12): \(X_i=B_iY_i\), particle mass-number/baryon fraction.
- Eqs. (4.20)-(4.21): \(A_{\rm av}\) and \(Z_{\rm av}\) are abundance-weighted means over a nuclear group.
- Eq. (4.22): \(Y_I=(\sum_j A_jY_j)/A_{\rm av}=\sum_jY_j\), the combined number fraction of the nuclear group.

Consequently, under standard CompOSE semantics:
\[
X_h=A_{\rm av}Y_I.
\]

The raw \(Y_I\) and the group baryon/mass fraction \(X_h\) are not interchangeable.

### 2.2 Native verifier semantics

Source: \`fortran_eos/src/compose_eos_verify.f90\`, C6-run identity blob SHA \`f1aa41c4a6e13ef7ea31c5f3c2329139b273b2d0\`.

The verifier evaluates:
\[
r_{C3}=Y_p+\frac{Z_{\rm av}}{A_{\rm av}}Y_h-Y_e-Y_\mu,
\]
\[
r_{C4}=Y_n+Y_p+Y_h-1.
\]

Default composition tolerance is \(2\times10^{-5}\).

These equations make the verifier's \(Y_h\) a baryon/mass fraction of the nuclear group, not its group number fraction. C3 is charge neutrality; C4 is baryon normalization. A proposed mapping must satisfy both equations on each row, not just preserve a field or fit one closure.

The reader maps CompOSE particle-pair identifiers as follows:
- 10: neutron
- 11: proton
- 0: electron
- 1: muon

The code's particle-pair mapping was inspected in \`fortran_eos/src/compose_eos_io.f90\`. Do not assume the nuclear-group quadruple's last field uses the same representation as verifier \`Yh\`.

## 3. Provider evidence and phase labels

### 3.1 SPG(M1)

Official [SPG(M1) data sheet](https://compose.obspm.fr/download/1D/Crust/M1/eos.pdf):
- Submitted by Luigi Scurto; data sheet dated 26 February 2024.
- Describes a BSk22 outer crust below \(n_B=2\times10^{-3}\,\mathrm{fm}^{-3}\), and an inner-crust EoS computed with a compressible liquid-drop approach.
- States that \`eos.compo\` contains four particle pairs (neutrons, protons, electrons, muons) and one heavy-nucleus quadruple.
- Explicitly labels \(Iphase=0\) outer crust and \(Iphase=7\) inner crust.

The data sheet does not, by itself, explain any non-standard normalization of the nuclear-group abundance field.

Original-work references listed by the data sheet include the BSk22 outer-crust paper, the compressible liquid-drop framework of Baym-Bethe-Pethick, the GDFM model paper, and an in-preparation Scurto/Pais/Gulminelli work. The agent must verify the exact original source/provenance and avoid treating an in-preparation citation as a fully accessible published specification.

### 3.2 GPPVA(FSU2)

Official [GPPVA(FSU2) data sheet](https://compose.obspm.fr/download/1D/Crust/FSU2/eos.pdf) and [CompOSE model entry](https://compose.obspm.fr/eos/209):
- The crust combines a BsK22 outer crust with an inner crust calculated in a self-consistent Thomas-Fermi approach using FSU2.
- \`eos.compo\` contains four particle pairs (neutrons, protons, electrons, muons) and one heavy-nucleus quadruple.
- The data sheet explicitly labels \(Iphase=0\) outer crust, \(Iphase=7\) inner crust and \(Iphase=6\) homogeneous matter.
- Its extra \`eos.thermo\` quantity describes pasta geometry: 0 outer crust, 1 droplets, 2 rods, 3 slabs, 4 tubes, 5 bubbles, 6 core.

Thus \(Iphase=7\) means inner crust in both selected provider data sheets. The fact that GPPVA uses a separate geometry code in \`eos.thermo\` must not be conflated with the \`Iphase\` code in \`eos.compo\`.

Original references identified by the provider include Grill et al., Phys. Rev. C 90, 045803 (2014); Chen and Piekarewicz, Phys. Rev. C 90, 044305 (2014); Providência et al., Front. Astron. Space Sci. 6, 13 (2019); and the BsK22 outer-crust work. The agent should retrieve and inspect the primary papers/data-generation descriptions where available.

## 4. Pinned raw source identity and representative rows

Source authority commit:
\`a09ec30d74b363652d2d1aefe990557c468679dd\`.

Raw \`eos.compo\` blob identities:
- SPG_M1: \`ab79da517271bea1f7de37294e87e019de611699\`
- GPPVA_FSU2: \`fa1ad2b809eb1c95825e2ee6eb2650cf3f248693\`

Manifest metadata:
- SPG_M1: 365 rows; \(Iphase=0\): 217; \(Iphase=7\): 148; boundary knots 217/218.
- GPPVA_FSU2: 269 rows; \(Iphase=0\): 217; \(Iphase=7\): 52; boundary knots 217/218.
- Both sources contain four particle pairs and one nuclear-group quadruple per row.
- GPPVA_FSU2 has \(Nadd=1\); its extra-column histogram is 0:217, 1:52. Do not treat that extra thermodynamic geometry code as the \`Iphase\` code.

### 4.1 Outer-crust representative row

Both sources begin with a row equivalent to:
\[
Iphase=0,\quad
Y_n=0,\quad Y_p=0,\quad Y_e=0.46429,\quad Y_\mu=0,
\]
\[
A_{\rm av}=56,\quad Z_{\rm av}=26,\quad Y_I=1.
\]

Current C5 internal mapping \(Y_h=Y_I=1\) yields:
\[
r_{C3}\simeq \frac{26}{56}-0.46429 \simeq 0,
\qquad
r_{C4}=1-1=0.
\]

The standard CompOSE product \(A_{\rm av}Y_I=56\) cannot be a normalized baryon/mass fraction for this row. This is a direct sign that the provider's stored raw value cannot be blindly used as standard \(Y_I\) with the manual's normal conversion, or that the source's composition representation requires further explanation.

### 4.2 First inner-crust row: SPG_M1

Raw row 218:
\[
Iphase=7,\quad Y_n=Y_p=Y_\mu=0,\quad
Y_e=0.06087,
\]
\[
A_{\rm av}=82.107,\quad Z_{\rm av}=24.387,\quad Y_I=0.06087.
\]

Current mapping \(Y_h=Y_I\) yields nonzero C3/C4 residuals. Independently of the stored \(Y_I\), C4 with \(Y_n=Y_p=0\) requires \(Y_h=1\), whereas C3 requires
\[
Y_h=\frac{A_{\rm av}}{Z_{\rm av}}Y_e\approx 0.205.
\]
These are mutually inconsistent. Thus, for this representative row, merely choosing a different scalar mapping of \(Y_I\) cannot satisfy both current closure equations simultaneously. The agent must test whether the raw pair/group fields, their interpretation, or the verifier's closure assumptions are incomplete or inapplicable for this provider row.

### 4.3 First inner-crust row: GPPVA_FSU2

Raw row 218:
\[
Iphase=7,\quad Y_n=Y_p=Y_\mu=0,\quad
Y_e=0.04865,
\]
\[
A_{\rm av}=740,\quad Z_{\rm av}=36,\quad Y_I=0.04865.
\]

C4 requires \(Y_h=1\). C3 independently requires
\[
Y_h=\frac{A_{\rm av}}{Z_{\rm av}}Y_e\approx 1.000.
\]
This row is therefore close to closure under \(Y_h=1\), but the stored raw \(Y_I=0.04865\) is not that value. This is an important contrast with SPG_M1 and is not grounds for adopting a source- or phase-specific correction without provider evidence and complete row-wise validation.

## 5. C6 native execution evidence

Evidence directory:
\`artifacts/sharp_maxwell/C6_NATIVE_EXECUTION_20261110T093714Z/\`.

Execution head:
\`43eea8bd9bee7ebb222e3c66babf5078b398d7ae\`.

Toolchain: gfortran 13.3.0, Python 3.12.3.

Verified:
- Official native build: PASS.
- Exact committed pytest file: 6/6 PASS.
- Native reader, writer, re-reader and verifier executed for both sources.
- Native read/write/read round-trip: PASS.
- Source hashes: 18 files unchanged.
- Legacy output hashes: 234 files unchanged.
- All 65 pair records unchanged and remain \`BLOCKED_BEFORE_ROOT_SEARCH\`.

Scientific composition checks:
- SPG_M1: C3/C4 fail on all 148 \(Iphase=7\) rows; C3 maximum absolute residual \(0.042790704812013594\); C4 maximum absolute residual \(0.9779100000000001\).
- GPPVA_FSU2: C3/C4 fail on all 52 \(Iphase=7\) rows; C3 maximum absolute residual \(0.04628324324324324\); C4 maximum absolute residual \(0.98856\).
- All \(Iphase=0\) rows pass C3/C4.

Interpretation: the native path and round-trip work, but physical composition semantics remain unresolved. The current \`Yh=YI\` mapping is not a scientific PASS.

## 6. Competing hypotheses to investigate

The local agent must test these alternatives rather than select one by intuition.

### H1. Provider-specific non-standard abundance convention

The providers may have populated the last nuclear-group field using a provider-specific quantity rather than standard CompOSE \(Y_I\), or may have used a convention that differs by source and/or phase. A provider-side export step, preprocessing script, or model Data Sheet appendix may document this.

### H2. Source composition data are inconsistent or incompletely populated

One or more composition fields may be wrong, omitted, or populated with an unintended quantity in some phases. Compare with original provider model outputs and published composition profiles. Check whether \(Y_n,Y_p,Y_e,Y_\mu,A_{\rm av},Z_{\rm av},Y_I\) independently correspond to the same grid point and same definition.

### H3. The selected native C3/C4 closure is not applicable to these provider outputs as encoded

The mathematical equations themselves are standard charge/baryon closure for a complete composition, but the model's one-group representation or omitted components could invalidate the current input interpretation. This must be established from the provider model and data schema, not assumed as a way to waive failing checks.

### H4. Parser/index/alignment defect

The C6 raw parser and native reader may agree on row ordering while still sharing an incorrect interpretation of a field or association. Verify the full record structure against the CompOSE file-format manual, grid keys \((i_T,i_{n_b},i_{Y_q})\), count fields, particle IDs, group IDs and ordering. Do not rely only on positional row equality.

## 7. C6B investigation DAG

\`\`\`mermaid
flowchart TD
    A["C6B-0: Freeze source/code identity"] --> B["C6B-1: Read CompOSE spec and both Data Sheets"]
    B --> C["C6B-2: Retrieve original provider papers and data provenance"]
    A --> D["C6B-3: Independently parse every raw eos.compo row"]
    C --> E["C6B-4: Reconstruct intended provider composition semantics"]
    D --> F["C6B-5: Compute candidate C3/C4 closures row by row"]
    E --> F
    F --> G["C6B-6: Compare source outputs with provider equations/data"]
    G --> H{"Unique, evidence-backed interpretation?"}
    H -->|Yes| I["C6B-7: Propose minimal code/data contract change"]
    H -->|No| J["STOP: PROVIDER_SEMANTICS_UNRESOLVED__BLOCKED"]
    I --> K["C6B-8: Independent review and explicit approval gate"]
\`\`\`

### Required tasks

1. Freeze the exact repository head, C5 code blob, verifier blob, raw source blob SHAs and C6 evidence identity. Current known branch head at task handoff is \`3f9531e6796a849b6f77b193bea2ecece74348f5\`; fetch it again before work because it may have advanced.
2. Read the CompOSE Reference Manual's definitions and record exact equation/page references. Do not conflate group number fraction \(Y_I\), group mass fraction \(X_h\), electron fraction \(Y_e\), or average \(A/Z\).
3. Inspect both official provider data sheets and retrieve the original provider papers, supplementary materials, repository data products and/or export scripts. The GPPVA(FSU2) record is associated with the Thomas-Fermi FSU2 model and the references above; SPG(M1) uses the BSk22 outer crust plus a compressible-liquid-drop inner crust. Trace the actual source of the exact submitted tables.
4. Independently parse every row of both raw \`eos.compo\` files directly from pinned source blobs. Validate count fields and keys; retain exact raw columns. Confirm that the raw reader's positional interpretation matches the official schema.
5. Build a per-row audit table with at least: source ID, row, \(Iphase\), any extra geometry code, particle IDs and fractions, \(A_{\rm av}\), \(Z_{\rm av}\), raw last quadruple field, \(Y_e+Y_\mu\), candidate \(Y_h\) values, C3 residual and C4 residual.
6. Evaluate at least these candidate interpretations without promoting any:
   - direct field \(Y_h=Y_I\);
   - standard CompOSE conversion \(Y_h=A_{\rm av}Y_I\);
   - values independently implied by C3 and C4;
   - any provider-specific formula only if a primary source documents it.
7. Compare all row-wise values with the provider's original model output or publication. Identify whether the first 217 outer-crust rows and the inner-crust rows use distinct production pipelines or definitions.
8. Decide separately for SPG_M1 and GPPVA_FSU2 whether:
   - the raw data conform to standard CompOSE;
   - the provider uses a documented non-standard representation;
   - the table is internally inconsistent or incomplete;
   - the current verifier's input contract is inappropriate.
9. If a unique interpretation is established, present the minimal proposed mapping/API contract and tests. Do not implement it in this investigation unless the user or project owner explicitly authorizes the repair after review.
10. Preserve all evidence and source identities. Use a machine-readable JSON summary and a human-readable report. Do not copy or modify the original raw EOS data.

## 8. Acceptance criteria

C6B is **resolved** only if all of the following are satisfied:

- Primary-source evidence documents the intended physical meaning of the stored nuclear-group field for each provider and relevant phase.
- The raw file format and every key/column association are independently verified.
- A complete row-wise recomputation of both C3 and C4 is available for both full source tables.
- The explanation accounts for both outer-crust and inner-crust behavior; no ad hoc phase switch is introduced solely to make closure pass.
- Proposed mapping(s) satisfy the provider's documented physics and reproduce the published/provider composition where a comparison is possible.
- Independent review can reproduce the conclusion from pinned source files and documented equations.

If these criteria are not met, report \`PROVIDER_SEMANTICS_UNRESOLVED__BLOCKED\`. Lack of documentation is not permission to guess.

## 9. Strict scope and stop conditions

Allowed:
- read-only research on the two selected source tables;
- inspection of official CompOSE documentation, original provider papers and data provenance;
- independent parsing and numerical closure recomputation;
- a proposed minimal remediation design.

Not allowed:
- modifying raw \`eos.nb\`, \`eos.thermo\`, \`eos.compo\`, or any source file;
- changing the current reader/writer mapping, verifier equations, or tolerances;
- inventing provider-specific phase-dependent conversions;
- changing any of the 65 pair records or their statuses;
- broadening the investigation to other EOS sources;
- executing a sharp-Maxwell root search, 65-pair sweep, TOV or Cowling;
- entering C7 whole-interval scientific verification before this blocker is resolved.

On any missing or contradictory provider evidence, stop with the exact unresolved point and the minimal evidence needed to decide it.

## 10. Terminal project state

Until C6B is resolved and the proposed semantics are explicitly accepted, preserve:

\`\`\`text
C6: C6_PARTIAL
C6B: OPEN / BLOCKED_REVIEW_REQUIRED
C7: NOT_RUN / NOT_RECOMMENDED
65 candidate pairs: BLOCKED_BEFORE_ROOT_SEARCH x65
legacy maxwell: LEGACY_NON_MAXWELL_READ_ONLY
sharp_transition_search: NOT_RUN
tov_cowling: NOT_RUN
\`\`\`

This review is a semantics/provenance investigation, not a mechanism for weakening scientific gates.
