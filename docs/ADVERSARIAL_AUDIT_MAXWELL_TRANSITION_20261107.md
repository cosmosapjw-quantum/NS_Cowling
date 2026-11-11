# Adversarial audit: native Fortran Maxwell transition v1

**Audit date:** 2026-11-07  
**Repository:** cosmosapjw-quantum/NS_Cowling  
**Target branch:** research/fortran-eos-gap-fill-v1  
**Pre-change code baseline:** 2eeae45c4c2e178b7577f674231ae400f3849806  
**Scope:** the newly added one-conserved-charge, locally neutral Maxwell construction and the code/documentation it touches. This audit does not certify the entire historic EOS generator or the stellar phase-interface solver.

## 1. Audit protocol and evidence rule

The review uses the requested sequence as an explicit, artifact-oriented checklist:

1. **Self-discover:** inspect the live branch tree, current Fortran EOS types, reader/writer, verifier, CLI, Python wrapper, parity matrix, tests and legacy capability declarations.
2. **Step-back:** separate well-posed work from underdetermined work. A single barotrope does not contain the independent charge densities or volume fraction needed for a Gibbs construction; a composition-free table cannot prove local electric neutrality.
3. **CRAG:** ground data conventions and physical assumptions in the CompOSE reference manual and independent phase-equilibrium literature, rather than treating legacy Python behaviour as authority.
4. **Metacognitive self-ask:** challenge the premise: are the supplied branches true locally neutral beta-equilibrated phases? Is the chemical potential shared? Are pressure, density and energy density derived consistently? Does the result preserve the source composition? Can a downstream TOV solver cross the interface correctly?
5. **CoVe:** independently re-derive and symbolically check the interpolation identities, sound-speed formula and latent-energy relation. Then connect each formula to a corresponding code gate or regression test.
6. **Adversarial self-ask:** attempt counterexamples involving multiple pressure crossings, tangential crossings, a missing neutrality assertion, acausal interpolation, non-overlapping chemical-potential domains, non-Euler-consistent input, unsupported provider-specific columns and serialization of a composition-free phase jump.
7. **CCoT-style claim-to-evidence trace:** map each implementation claim to its equations, source files, tests and final CI logs. A code path without a passing exact-HEAD build/test run is not labelled verified.
8. **PDR-style decision record:** preserve provenance, decisions, changes, blockers and the bounded stop condition. No claim of a native subagent review is made; native subagent dispatch was not available in this execution environment.

The terms above identify review activities and their durable outputs. They do not substitute for a reproducible derivation, test log, source reference or code diff.

## 2. Findings before changes

| Severity | Finding | Why it mattered | Resolution |
|---|---|---|---|
| Critical | The historical maxwell smooth-matching token could be confused with a physical Maxwell construction. | The old routine has no density jump or latent heat. | A separate maxwell-transition kernel and API were introduced; the old token remains explicitly documented as a smooth compatibility bridge. |
| Critical | Local charge neutrality was not encoded in composition-free input data. | A Maxwell construction assumes local neutrality in each phase; the assumption cannot be proved from a one-dimensional barotrope. | The caller must explicitly pass assume_local_charge_neutrality=True; the report records caller_confirmed. |
| High | A Gibbs mixed phase was not determined by the available data (Constantinou et al. 2023). | Global charge neutrality needs independent charge densities, a charge chemical potential and a volume fraction. | Gibbs is not implemented and is explicitly excluded. |
| High | Dense point sampling alone is not a proof of mechanical stability or causality between knots. | A cubic interpolant can violate an inequality between sample points. | The Hermite polynomial is analysed algebraically over each interval: n(mu) and mu n'(mu)-n(mu) are quadratic; n'(mu) is affine. |
| High | A fixed-grid pressure-root scan can miss a narrow pair of crossings or a tangential root. | The pressure difference is cubic on each union-knot interval. | Stationary points from its quadratic derivative are isolated first; each monotone subinterval is checked and bracketed. |
| High | Provider-specific Nadd fields had no generic interpolation law. | Discarding or interpolating them would invent model semantics. | V1 rejects nonzero Nadd. |
| High | The phase-boundary derivative mask depended on composition availability. | A phase discontinuity exists even when eos.compo is absent. | C7-C9 derivative stencils now mask both rows adjacent to every Iphase change. |
| High | A phase-labelled barotrope without phase composition could be exported as if it were a complete CompOSE table. | The output does not contain phase-specific charge fractions or eos.compo. | Both the Python wrapper and native Fortran writer refuse this export. |
| Medium | Some negative API tests did not pass the required neutrality assertion. | They failed at the API precondition instead of exercising the intended native validation. | Tests were corrected and a separate missing-assertion test was added. |
| Medium | The Python writer guard displaced the function docstring. | This broke introspection/documentation behaviour. | The docstring was moved back to the first statement. |

## 3. Five bounded review-and-fix cycles

### Cycle 1: physical scope and conventions

**Questions:** What exactly is conserved? Are both phases locally neutral? Is beta equilibrium represented? Is the historic maxwell method actually a first-order transition?

**Evidence:** the CompOSE Reference Manual v3.01 (Typel et al. 2022); Hempel et al. (2009); Bhattacharyya et al. (2010); the existing six-method smooth matching code and its report fields.

**Change:** Introduced a separate Maxwell routine for one globally conserved charge and local neutrality in each supplied phase. The local-neutrality assumption is explicit. Gibbs and composition reconstruction were excluded rather than inferred from absent columns.

**Status:** Design and source-level scope fixed. External physical assumptions are literature-supported; local neutrality remains caller-supplied provenance.

### Cycle 2: thermodynamic closure and causal representation

**Questions:** Does the interpolation preserve dP=n_B dmu_B? Does the energy density obey epsilon=mu_B n_B-P? Is 0<c_s^2/c^2<=1 enforced throughout each interval?

**Derivation:** With n_B=P'(mu_B) and epsilon=mu_B n_B-P,

\[
\epsilon'=\mu_Bn'_B,\qquad
\frac{c_s^2}{c^2}=\frac{n_B}{\mu_Bn'_B}.
\]

Thus positive n_B, positive n'_B, and
mu_B n'_B-n_B>=0 are the required conditions.

**Change:** Added the cubic-Hermite P(mu_B) representation and analytic polynomial-extremum checks. The 33-point diagnostic is no longer the acceptance proof.

**Independent check:** Wolfram Language symbolic evaluation returned zero for P'-n, epsilon'-mu n', and P'/epsilon'-n/(mu n') for a generic cubic P(mu) with n=P'.

**Status:** Algebraically derived and symbolically checked. Native build and analytic regression verified on the exact Maxwell feature revision.

### Cycle 3: pressure-root topology and transition identity

**Questions:** Can the search miss an even number of crossings? Can a tangent crossing be found? Are multiple candidate crossings selected silently? Does the density jump satisfy the latent-energy identity?

**Change:** The union of both knot grids partitions the chemical-potential domain into intervals where each phase has one Hermite polynomial. On each interval delta-P is cubic and delta-n=delta-P' is quadratic. The derivative roots partition it into monotone subintervals, on which sign changes and endpoint roots are tested. Exactly one unique coexistence root is required. The quadratic formula was written in a cancellation-resistant form.

**Output check:**

\[
P_L(\mu_t)=P_H(\mu_t),\quad
\mu_L=\mu_H=\mu_t,\quad
\Delta\epsilon=\mu_t\Delta n_B.
\]

**Status:** Mathematical construction implemented and the analytic EOS regression fixture passed exact-head native CI.

### Cycle 4: phase boundary, metadata and export safety

**Questions:** Does the verifier differentiate across the density jump? Can phase labels survive without composition? Could a user write an incomplete phase EOS as a complete CompOSE product?

**Change:** The verifier masks derivatives adjacent to any phase-label change. The native writer and Python wrapper reject a phase-discontinuous table if the phase composition is absent. The output labels the two branches but does not claim to reconstruct eos.compo.

**Status:** Source-level safeguards and regression tests added. CI is the gate for implementation verification.

### Cycle 5: API, documentation, tests and release gate

**Questions:** Do negative tests reach their intended failure conditions? Is the assumption visible in both API and CLI? Does the documentation state what is not implemented? Does the source archive include the exact tested tree?

**Change:** Corrected API test flags; added explicit-neutrality rejection, composition-rejection, non-overlap, acausal-input, Q3-Euler mismatch, Python export refusal and native CLI export refusal tests. Added a self-contained derivation/implementation note and an updated conservative feature matrix.

**Status:** The first exact-HEAD CI run on `6ffdd62080dd8313e6c890cf6a6261e68f5c864d` compiled the complete Fortran executable successfully. The suite then reported **53 passed, 1 failed**: the analytic Maxwell test expected `low_phase_rows`, but the CLI JSON omitted that field although the native report type contained it. This is a Python/native report-schema serialization defect, not a failure of the coexistence equations. The CLI JSON has been extended with `low_phase_rows`, `high_phase_rows`, and phase-resolved sound-speed extrema. The Maxwell feature passed both push and PR CI on commit 37f7fdcba998b05f73f8ab2ee94e85c24ea37b19: 54 tests passed, the then-current 62-entry manifest passed, and source archive creation/upload succeeded.

## 4. Claim-to-evidence table

| Claim | Evidence location | Required verification |
|---|---|---|
| CompOSE Qi meanings and units | CompOSE Reference Manual v3.01, arXiv:2203.03209 | External reference |
| Maxwell versus Gibbs assumptions | Hempel et al. (2009), PRD 80, 125014 | External reference |
| dP/dmu_B=n_B and epsilon=mu_B n_B-P | compose_eos_phase_transition.f90 | Code review + analytic EOS test |
| Full-interval n>0, n'>0, causal bound | quadratic_interval_minimum and validate_phase in same module | Native execution |
| Exactly one pressure-equality root | find_coexistence | Root-topology tests and native execution |
| Density jump and latent-energy identity | construct_maxwell_transition, analytic test | Native execution |
| Phase-boundary finite-difference masking | compose_eos_verify.f90 | Existing verifier tests + CI |
| No incomplete CompOSE export | compose_eos_io.f90, Python write_compose | Direct CLI and Python tests |
| No Gibbs or TOV interface claim | parity matrix and this note | Documentation consistency |

## 5. Literature and source provenance

- S. Typel, M. Oertel and T. Klähn (2022), *CompOSE Reference Manual v3.01*, [arXiv:2203.03209](https://arxiv.org/abs/2203.03209). Primary convention source for Q1,...,Q7, the lepton flag and separation of thermodynamic and composition tables.
- M. Hempel, G. Pagliara and J. Schaffner-Bielich, *Conditions for phase equilibrium in supernovae, protoneutron, and neutron stars*, Physical Review D, vol. **80**, article 125014, 2009, [doi:10.1103/PhysRevD.80.125014](https://doi.org/10.1103/PhysRevD.80.125014), [arXiv:0907.2680](https://arxiv.org/abs/0907.2680). Equilibrium conditions with globally conserved charges and locally constrained charge fractions.
- A. Bhattacharyya, I. N. Mishustin and W. Greiner, *Deconfinement phase transition in compact stars: Maxwell vs. Gibbs construction of the mixed phase*, Journal of Physics G, vol. **37**, article 025201, 2010, [doi:10.1088/0954-3899/37/2/025201](https://doi.org/10.1088/0954-3899/37/2/025201).
- *Framework for phase transitions between the Maxwell and Gibbs constructions*, Physical Review D, vol. **107**, article 074013, 2023, [doi:10.1103/PhysRevD.107.074013](https://doi.org/10.1103/PhysRevD.107.074013).

## 6. Blockers and explicit non-claims

1. **Initial CI defect preserved:** commit 6ffdd620 built but failed one test (53 passed, 1 failed) because the CLI JSON omitted low_phase_rows and high_phase_rows. Commit 37f7fdc added these keys and passed the full 54-test, checksum and packaging gates.
2. **No local compiler run:** the connected Desktop Commander endpoint reported its monthly usage limit. This was not retried. GitHub Actions is the available compiler/test environment.
3. **No native subagent dispatch:** no tool for independent native subagent execution was available. The review cycles are direct source/math/test audits, not simulated subagent reports.
4. **Gibbs, composition reconstruction, full reconstruct_master, literature-specific EOS solvers, Seidov threshold and phase-interface TOV stability remain unverified or unported.**
5. **No merge to main:** the intended change is confined to the existing feature branch and its open PR.
6. **Initial regression outcome retained:** the first report-schema failure is preserved as a historical defect; it was corrected in the later exact-head commit and is not represented as a passing run.

## 7. Stop condition

The bounded v1 task is complete only when the exact feature-branch commit passes:
- Fortran 2018 compilation with the repository CI flags;
- the complete targeted frontend/direct-provider test suite, including the analytic Maxwell fixture and writer-refusal tests;
- both 58-plus tracked SHA-256 manifest entries after the update;
- source archive generation/upload.

A failed gate returns to the relevant bounded cycle. A passing gate establishes implementation verification for this v1 contract only; it does not close the separately listed physics backlog.


## 8. Writing and citation-tool check

The two technical documents were checked with the Academic Writing Toolkit. Both passed its conservative British-English spelling check with zero findings. Its paragraph-logic pass flagged short standalone equation lead-ins and concise exclusion statements; these were retained where they separate displayed mathematics, code interfaces or explicit non-claims. The citation audit did not yield a usable PASS: it reported `notes-source-malformed` for a Harvard `**Source**` line returned by the toolkit's own reading-note template generator, followed by phantom-citation warnings. We therefore treat that automated citation result as **INCONCLUSIVE**, not as a clean audit. Bibliographic entries and in-text author-year citations were instead cross-checked manually against the cited CompOSE manual and DOI/arXiv records. This tool-schema inconsistency is retained here rather than concealed.


## 9. Exact Maxwell release evidence

On commit 37f7fdcba998b05f73f8ab2ee94e85c24ea37b19, push CI [run 42381968318](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42381968318) and PR CI [run 42381971708](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42381971708) both completed successfully. The push workflow reported 54 passed in 21.41 s, 62/62 manifest entries, source archive packaging and artifact upload. Later master-reconstruction commits continued to build and test on the same feature branch; see the separate master reconstruction audit.


## 9. Exact Maxwell release evidence

On commit 37f7fdcba998b05f73f8ab2ee94e85c24ea37b19, push CI [run 42381968318](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42381968318) and PR CI [run 42381971708](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42381971708) both completed successfully. The push workflow reported 54 passed in 21.41 s, 62/62 manifest entries, source archive packaging and artifact upload. Later master-reconstruction revisions build and test on the same feature branch; see the separate master reconstruction audit.


## 10. Cycle 6: CompOSE Il=0 Q3 convention (2026-11-11)

**Claim under audit:** the true one-charge Maxwell kernel accepts only source phases whose stored Q3 agrees with the Euler baryon chemical potential after the CompOSE lepton-flag convention is applied, and writes all retained/output rows in that same convention.

**Source-level defect:** validate_phase formed mu_source=mn*(1+Q3) and compared it directly to mu_E=(epsilon+P)/n_B for both Il=0 and Il=1. The module permits either flag. The repository C2 verifier and reconstruct_master instead apply the Il=0 correction mu_Q3=mu_E-Yq*mn*Q4. The set_transition_row and copy_source_row output paths also wrote Q3 as if Il=1 for both. A valid Il=0 phase could therefore be rejected, and an Il=0 output would fail ordinary C2 verification.

**Required convention:**

\[
\mu_E=\frac{\epsilon+P}{n_B},\qquad
m_n(1+Q_3)=
\begin{cases}
\mu_E-Y_qm_nQ_4,&I_l=0,\\
\mu_E,&I_l=1,\ Q_5\simeq0.
\end{cases}
\]

The coexistence coordinate remains Euler mu_E for both branches. This correction does not create another conserved charge, does not change Maxwell equilibrium conditions, and does not authorize Gibbs/global-neutrality construction.

**TDD evidence:** the new analytic two-phase Il=0 regression was run before the native edit and failed at the intended input gate (provider Q3 chemical potential disagrees with the Euler value). The minimal native fix adds Yq*mn*Q4 to stored Q3 before source/Euler comparison, and subtracts Yq*Q4 from Q3 when writing both coexistence rows and copied source rows. After the fix, the targeted regression passed and check_consistency returned C2 pass and overall pass.

**Verification:** GNU Fortran 14.2.0 build passed; native polynomial-extrema test passed; the complete documented frontend/pipeline suite collected and passed 68 tests in 13.57 s. The exact source bytes were executed from the PR #1 GitHub Actions source artifact exported into a local directory without .git metadata. This is an artifact-based local execution, not a committed Git worktree replay. The repository SHA256 manifest was refreshed after the final tracked edits and checked separately before publication.

**Scope limits retained:** only the two supported cold branches Il=0 and Il=1,Q5≈0 are covered. Nonzero-Q5/lepton-rich states remain rejected; composition-bearing/Nadd inputs remain outside Maxwell v1; no microscopic charge neutrality is inferred; no Gibbs mixed phase, Seidov verdict, phase-interface TOV stability or Cowling spectrum is established. The original single-phase reconstruction and six smooth crust-core methods were not altered.
