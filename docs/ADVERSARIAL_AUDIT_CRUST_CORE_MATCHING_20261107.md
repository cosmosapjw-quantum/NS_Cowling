# Adversarial audit: native smooth crust-core EOS matching

**Repository:** cosmosapjw-quantum/NS_Cowling  
**Target branch:** research/fortran-eos-gap-fill-v1  
**Baseline HEAD:** 6b4cabc05e60bbe159831eec2397d45df66aa5bb  
**Scope:** native port of six historical smooth thermodynamic matching methods; not full EOS-generator parity.  
**Review budget:** maximum five bounded review/fix cycles.  
**Main branch policy:** no merge to main.

## Evidence vocabulary

This audit distinguishes source inspection, literature context, mathematical derivation, symbolic verification, synthetic numerical testing, native differential tests and CI evidence. Repeated self-review is not an independent reviewer. No subagent review is claimed unless a native agent execution is available and recorded.

## Stage mapping

The requested procedure is treated as an inspectable evidence workflow. Each stage below is tied to a repository source, mathematical check, adversarial test or durable decision artifact:

- **Self-discover:** identify the active implementation, tests, feature matrix, legacy reference and branch state.
- **Step back:** define the physical problem narrowly and reject scope inflation into full generator parity.
- **CRAG:** retrieve the legacy source, stored benchmark metadata, official CompOSE documentation and EOS interpolation literature.
- **Metacognitive self-ask:** search for silent information loss, invalid input states, ambiguous transition labels and criteria that may pass while the output remains unphysical.
- **CoVe:** independently derive the thermodynamic formulae and check the algebra symbolically.
- **Adversarial self-ask:** search for spline-index, endpoint, fallback, root-bracketing, monotonicity and false-Maxwell errors.
- **CCoT:** compare analytically distinct pressure--chemical-potential, free-energy, CSS and pressure-plateau formulations.
- **PDR:** record the decision, evidence, limitations and stop condition.


## Bounded review/fix iterations

### Iteration 1: repository and scope discovery

- Re-read PR #1, the exact branch, native types/reader/writer/CLI/build script, frontend API, feature matrix and legacy matching source.
- Confirmed six historical tokens and a strict thermo-only mode.
- Confirmed the historical maxwell token constructs a smooth near-pressure plateau without density jump or latent heat; it is not a physical Maxwell coexistence construction.
- Fixed v1 scope to strict thermo-only. Composition reconstruction and literature_beta_eq are excluded.
- Native subagent review: **NOT_RUN**; no native agent-delegation capability was available.

### Iteration 2: thermodynamic derivation and source grounding

For cold one-component matter,
\[
 \mu_B=F+Q_1,\quad dF=Q_1d\ln n_B,\quad dP=n_Bd\mu_B,
 \quad c_s^2=\frac{dP/dn_B}{d\epsilon/dn_B}.
\]
The pressure bridge therefore integrates dP=n_B dmu_B, whereas the free-energy bridge integrates dF=Q1 dln n_B. For the latter,
\[
 c_s^2=\frac{Q_1+dQ_1/d\ln n_B}{F+Q_1}.
\]
The quintic Hermite free-energy representation was differentiated independently with Wolfram Language; the difference between the analytic second derivative and the implemented coefficient expression simplified exactly to zero.

The legacy source was inspected for transition-window screening, one-sided PCHIP endpoint sound speeds, gauge alignment, all six bridges and output composition behavior. Official CompOSE documentation and relevant interpolation literature were retrieved. References and non-equivalence limits are in the self-contained implementation note.

### Iteration 3: implementation adversarial review and correction

The initial implementation review found a piecewise spline evaluator defect: after locating the active segment, the integral accumulator could use the wrong interval. It was corrected to preserve the active segment and sum exact primitives over all preceding segments before evaluating the current segment. This is important because the integral enforces the first law.

The cubic admissibility threshold was corrected to agree with the legacy canonical free-energy branch. Internal spline helpers were made private so the supported API is the matching routine, not its implementation detail.

### Iteration 4: synthetic numerical and negative tests

A synthetic cold barotrope exercised all six methods. The scratch native results showed increasing density and pressure, positive energy density, bridge sound speed inside the stable causal interval, and endpoint First-Law residuals of approximately 1.3e-13 for the near-plateau bridge, 8.3e-14 for CSS and at or below round-off for the other four methods.

The fallback paths were separately forced in a local harness before internal helpers were made private. The cubic one-interior-knot profile met the requested integral. The adaptive C2 two-interior-knot profile preserved endpoint derivatives and met its integral. A subsequent attempt to relink that old harness after the helper became private failed at link time because the harness called an intentionally private symbol. This is an obsolete harness/API mismatch, not a production-code test failure; the supported six-method public routine remained the synthetic test target.

Input-state guards were added for T=0, beta-equilibrium metadata, zero entropy, Q6 approximately equal to Q7, and explicit consent before composition is discarded.

### Iteration 5: legacy differential tests, release checks and packaging

The latest pre-integration CI run (42374515991) on concurrent HEAD 0b27bed3aab8a6e39015c5984be6c3cabcc90f75 failed during compilation, before tests: the phase-local causal-stencil condition lacked a closing parenthesis, and free_per_baryon was declared scalar but assigned an array expression. Both compile defects were corrected without reverting the phase-local crossing logic. The matching module, CLI dispatch, build order, Python wrapper and tests were then integrated on top of that exact HEAD.

The first integrated CI attempt (42376763661) compiled the complete native executable successfully but ended with 34 passed and 9 failed. Six failures came from the legacy package's eager plotting import requiring matplotlib, so the differential harness now loads only the legacy EOS I/O and matching modules without executing the package-level plotting imports. One failure exposed the prior no-argument NativeFeatureUnavailable API contract, which is now explicitly preserved when one of the two EOS inputs is missing. One test failed only because the Python error wording used an underscore rather than the documented hyphenated phrase; that message is normalized. The remaining failure exposed a real metadata parser defect: Fortran L1 logicals are written as T/F, but Python parsed only .true./.false., making the string F truthy. The parser now handles both encodings and has a focused regression test. The six numerical differential comparisons subsequently ran. Five methods passed native-versus-legacy array comparisons; the near-pressure-plateau `maxwell` compatibility token differed in bridge Q1 by up to 5.27e-5 relative. Inspection isolated the mismatch to quadrature policy: the legacy code uses 64-point Gauss--Legendre for the full plateau root and 16-point Gauss--Legendre for partial free-energy integrals, while the first native implementation used Simpson quadrature. Native v1 now reproduces the same two Gauss--Legendre orders. The standalone T/F metadata parser test also initially contained literal backslash-n characters instead of newline escapes; the fixture was corrected. The GL64/GL16 correction and the test-fixture newline repair were committed, then the final release gate passed on exact HEAD `3486e3ca9407ddc1362208b3d36f0c3154d78e6e`: native gfortran build PASS, **44 tests passed**, all 58 entries in `SHA256SUMS_v4.txt` PASS, and source-snapshot artifact upload PASS. Both push CI run 42377746709 and pull-request CI run 42377750637 completed successfully on that same HEAD. The audit also records the manifest helper defect and its correction: exact UTF-8 bytes were hashed with Wolfram SHA-256 over ByteArray and independently accepted by CI `sha256sum -c`.

**Bounded feature decision: IMPLEMENTATION-VERIFIED.** The six smooth strict-thermo-only matching methods are compiled, wrapped, differentially tested against the historic Python reference, and packaged. This does not confer complete EOS-generator parity: composition-aware matching, `literature_beta_eq`, full master-function smoothing/candidate selection and the other separately listed generator/TOV capabilities remain unported. The push workflow produced artifact `NS_Cowling-3486e3ca9407ddc1362208b3d36f0c3154d78e6e` (198,966,313 bytes; SHA-256 `97e3f62e03b8bab2a7d11cebb56b671fefe0f624c5791349fc9a4f630579fc5c`; expires 2026-12-07).

## Principal adversarial questions

| Question | Safeguard or disposition |
|---|---|
| Can a finite sentinel masquerade as a valid derivative? | Missing derivatives are not encoded with finite sentinels. |
| Can endpoint tests pass while the bridge is mechanically unstable inside? | A dense-grid check tests positive pressure/energy, monotonicity and 0<c_s^2<1. |
| Can the integral pass while the wrong spline segment is evaluated? | Segment selection was corrected; the integral residual is reported. |
| Can composition disappear silently? | Python and Fortran both require explicit drop_composition=True; output metadata records omission. |
| Can the method name overclaim a phase transition? | maxwell is recorded as a smooth near-plateau compatibility bridge, with zero density jump and Seidov inapplicable. |
| Can a selected window be mistaken for a convexity theorem? | Secant inequalities are documented as screening heuristics; the constructed bridge is checked separately. |
| Can a local bridge pass certify the whole EOS? | No; the full C1--C11 report is returned and can set diagnostic_only=True. |
| Can native v1 be mistaken for full Python parity? | The parity matrix labels the feature PARTIAL and lists composition and microscopic-model limitations. |

## Stop condition and decision

The bounded strict-thermo smooth matching increment satisfies its stop condition: the exact tested HEAD compiles, all 44 targeted frontend/direct-provider tests pass, the 58-entry SHA-256 manifest passes, and the exact source commit was packaged. It is implementation-verified within the stated scope. It is **not complete legacy EOS-generator parity**; composition-aware reconstruction, literature_beta_eq and the other distinct legacy capabilities still require separate implementations and validation.

Main remains unchanged. PR #1 remains open and unmerged pending an explicit user decision.

## Final integrated release evidence

GitHub Actions run 42377336346 on commit 4b7f6a983e3a7cbcd920dad7cdfee3caee7125dc completed the native gfortran build and all **44 tests passed in 24.51 s**. This includes the six native-versus-legacy matching-method differential tests, no-argument API compatibility, Fortran L1 metadata parsing, existing causal-completion regressions, and direct-provider/frontend tests. The checksum stage initially exposed a stale density-gap test digest and then the audit document's digest. The first attempted repair used an auxiliary JavaScript SHA-256 routine with an incomplete compression-constant table; that helper was discarded. Both exact UTF-8 byte streams were re-hashed with Wolfram SHA-256 over ByteArray, the two manifest copies were updated, and the final exact-HEAD CI independently accepted every manifest entry.


## Manifest hash calculation provenance

The first repair attempt used an auxiliary JavaScript SHA-256 routine whose compression-constant table was incomplete for long inputs. This did not affect source code or test execution, but it produced invalid expected digests for two long files. Those outputs were discarded. The exact UTF-8 byte streams were re-hashed with Wolfram's SHA-256 over ByteArray values and will be checked independently by the final CI `sha256sum -c` step. This audit records the helper defect rather than treating the mismatch as a repository or runtime failure.
