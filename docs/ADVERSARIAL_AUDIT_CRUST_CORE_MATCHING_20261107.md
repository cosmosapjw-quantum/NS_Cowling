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

The first integrated CI attempt (42376763661) compiled the complete native executable successfully but ended with 34 passed and 9 failed. Six failures came from the legacy package's eager plotting import requiring matplotlib, so the differential harness now loads only the legacy EOS I/O and matching modules without executing the package-level plotting imports. One failure exposed the prior no-argument NativeFeatureUnavailable API contract, which is now explicitly preserved when one of the two EOS inputs is missing. One test failed only because the Python error wording used an underscore rather than the documented hyphenated phrase; that message is normalized. The remaining failure exposed a real metadata parser defect: Fortran L1 logicals are written as T/F, but Python parsed only .true./.false., making the string F truthy. The parser now handles both encodings and has a focused regression test. The six numerical differential comparisons subsequently ran. Five methods passed native-versus-legacy array comparisons; the near-pressure-plateau `maxwell` compatibility token differed in bridge Q1 by up to 5.27e-5 relative. Inspection isolated the mismatch to quadrature policy: the legacy code uses 64-point Gauss--Legendre for the full plateau root and 16-point Gauss--Legendre for partial free-energy integrals, while the first native implementation used Simpson quadrature. Native v1 now reproduces the same two Gauss--Legendre orders. The standalone T/F metadata parser test also initially contained literal backslash-n characters instead of newline escapes; the fixture was corrected. The next exact-head CI run is the release gate for these two fixes.

**Matching-specific release gate remains pending.** The targeted regression suite compares all six methods against the legacy Python package using the repository GPPVA_FSU2H crust and GM_GM1 core. It checks transition indices, n_B, Q1, Q3, Q4, Q6 and Q7, source-knot retention, composition flags, monotonicity, bridge causality and the First-Law residual. The exact integrated commit must pass compilation, the differential suite and SHA-256 manifest validation before the feature is called implementation-verified. Prior density-gap or causal-completion CI is not evidence for this matching increment. The workflow also packages that exact source commit as a downloadable artifact.

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

The implementation may be reported as **native strict-thermo smooth matching implemented** only after the exact branch HEAD has passing compilation, differential tests and manifest verification. It may be called **complete legacy parity** only after composition-aware reconstruction, literature_beta_eq and all other distinct legacy capabilities are ported and validated. That larger condition is not met here.

Main remains unchanged. PR #1 remains open and unmerged pending an explicit user decision.

## Final integrated evidence before manifest repair

GitHub Actions run 42377336346 on commit 4b7f6a983e3a7cbcd920dad7cdfee3caee7125dc completed the native gfortran build and all **44 tests passed in 24.51 s**. This includes the six native-versus-legacy matching-method differential tests, no-argument API compatibility, Fortran L1 metadata parsing, existing causal-completion regressions, and direct-provider/frontend tests. The tracked manifest check then found one stale SHA-256 entry for `python_frontend/tests/test_native_density_gap.py`; the file was unchanged by this matching increment, but its manifest entry no longer matched the branch's actual bytes. The exact file SHA-256 was recalculated as `4a9f4bb6542a5a7e0b3fef4f6079b74578e33208d9d4ddc2c83de035ee891e12` and both manifest copies were updated. A final exact-HEAD CI run must still recheck the manifest and create the packaged source artifact.
