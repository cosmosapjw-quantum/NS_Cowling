# Adversarial audit record: Fortran EOS generator port

**Audit target:** `cosmosapjw-quantum/NS_Cowling`, feature branch `research/fortran-eos-gap-fill-v1`  
**Baseline before this continuation:** `6b4cabc05e60bbe159831eec2397d45df66aa5bb`  
**Audit contract:** no merge to `main`; no unverified build or test claims; distinguish source review, symbolic derivation, numerical test, and CI evidence.  
**Review limit:** at most five bounded review/fix iterations. A review does not recursively audit another review.

## 1. Scope decision

The repository already contained a native density-gap fill and still lacked several large EOS-generation capabilities. This continuation selects the native port of the legacy constant-sound-speed causal completion as one bounded increment. It does not claim whole-generator parity. Full master-function smoothing, six crust-core matching procedures, microscopic literature solvers, genuine Maxwell/Gibbs transitions, and a fully native causal profile remain deferred.

## 2. Audit vocabulary and evidence policy

The requested named adversarial sequence is treated as a checklist for externally inspectable work products: discover the actual code surface; step back to fix the scope and scientific question; retrieve and compare authoritative source material; ask what could invalidate the claim; verify equations independently; conduct adversarial implementation checks; compare the construction with competing representations; then record a decision and reproducible release evidence. This record reports findings and artifacts, not private chain-of-thought.

A claim is classified as one of: source-inspected, literature-supported, derived, symbolic-CAS-checked, numerically-tested, implementation-verified, CI-verified, unresolved, or blocked. Agreement between repeated reviews is not treated as evidence by itself.

## 3. Review/fix iterations

### Iteration 1: discover and step back

- Read the current PR, branch head, full repository tree, native Fortran EOS types/CLI/build script, Python EOS frontend, legacy `causality.py`, EOS feature-parity matrix, and existing targeted tests.
- Confirmed that legacy Python `causal_completion` is a separate high-density CSS continuation, while the six crust-core matching methods are a much larger and physically distinct feature.
- Fixed scope to a native implementation of the bounded CSS completion. Other parity gaps remain explicitly listed.

**Evidence:** repository source files at the baseline branch head; no local source checkout was available in the connected execution environment.

### Iteration 2: source retrieval and independent thermodynamic derivation

- Retrieved the legacy algorithm and relevant literature records for CSS EOS parameterizations and causal EOS representations.
- Derived the continuation from (dP/dn=s(\epsilon+P)/n) and (d\epsilon/dn=(\epsilon+P)/n).
- Independently differentiated the closed-form (P(n)) and \(\epsilon(n)\) expressions with Wolfram Language. The pressure and energy derivatives have ratio (s), and (mu_B=(\epsilon+P)/n) scales as (n^s).
- Distinguished the implementation from sharp-transition CSS constructions and from a unique microscopic prediction.

**Evidence:** Wolfram Language symbolic result; Alford et al. (2015), Lindblom (2018), and Fortin et al. (2016), listed in the implementation note.

### Iteration 3: implementation-level adversarial checks

Issues identified and addressed before release candidate:

- A finite sentinel such as `-huge()` is still IEEE-finite. It cannot be used to encode a missing derivative while later relying on `ieee_is_finite`. The native routine now carries a separate logical validity mask for knot sound speeds.
- The crossing interval must not use a PCHIP derivative stencil that crosses an `Iphase` boundary. The native routine checks its local support stencil and refuses a phase-edge-only crossing.
- A negative source (c_s^2) knot below the first causal-limit crossing is a mechanical-instability blocker, not an acceptable prelude to the CSS tail.
- The chemical-potential match is independently checked against the Euler value with the legacy (5\times10^{-3}) relative threshold.
- The generated output does not preserve microscopic composition or phase labels. This is explicit in both the output metadata and the documentation.

**Evidence:** source diff and new synthetic regression tests; CI execution pending until the final branch commit.

### Iteration 4: numerical and runtime verification

**Pending until the final branch commit's GitHub Actions run.** Record here the actual build, exact test count, artifact/checksum result, workflow run ID, and HEAD SHA. Do not substitute the previous density-gap CI run for evidence about this increment.

### Iteration 5: final package and independent documentation audit

**Pending until Iteration 4 passes.** Re-read the final source/docs at the exact HEAD, check the release package contents and manifest, review scientific claims against equations and citations, and record unresolved limitations. No release-complete label is admissible before this row is filled.

## 4. Current blockers and declared limitations

- The connected Desktop Commander device was online but its monthly execution quota was exhausted; no retry or reconnect was attempted. The container could not resolve GitHub for a local clone. The source changes are therefore being authored through GitHub Git objects and verified by GitHub Actions, not by claiming a local build.
- Native subagent tools were not yet used. If no host delegation capability is available, the final record must say native subagent execution was NOT_RUN; serial self-review is not to be relabelled as independent agent review.
- The native causality profile is not part of this bounded port.
- The remaining Python EOS generation and matching pipeline is not fully ported.
- A passing compiler/test run does not establish observational validity of the CSS choice. The tail remains a user-selected phenomenological model assumption.
