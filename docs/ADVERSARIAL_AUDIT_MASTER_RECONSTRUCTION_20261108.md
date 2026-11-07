# Adversarial audit: native Fortran master-function reconstruction v1

**Audit date:** 2026-11-08  
**Repository:** cosmosapjw-quantum/NS_Cowling  
**Branch:** research/fortran-eos-gap-fill-v1  
**Pre-change head:** 37f7fdcba998b05f73f8ab2ee94e85c24ea37b19  
**Validated code head:** be648addc3fc20a937051b5b623009c9ef2d7eb8  
**Scope:** native cold master-function reconstruction, composition-aware electron/muon closure, remnant-nucleus repair, outer-crust electron/lattice correction, Python/CLI integration and associated documentation. This audit does not certify a full EOS-generator port, Gibbs construction, Seidov stability, or phase-interface TOV/Love integration.

## 1. Review protocol and evidence standard

The requested adversarial sequence was converted into bounded activities with inspectable outputs.

1. **Self-discover:** inspected the live branch tree, CompOSE EOS type and metadata model, native reader/writer, existing gap-fill and matching routines, the legacy master reconstruction, consistency checker, frontend wrappers, tests, workflow and release manifest.
2. **Step-back:** separated master-function reconstruction from complete EOS generation. A supplied free-energy master function can define cold barotropic thermodynamics, but it does not determine equilibrium nuclear species, microscopic couplings, global smoothing policy, phase-transition surface tension or stellar radial stability.
3. **CRAG/source grounding:** used the CompOSE Reference Manual for the Q1-Q7 convention and the BPS/outer-crust literature for relativistic electrons and Coulomb-lattice terms. Repository-specific smoothing and dissolution rules are explicitly treated as numerical implementation choices, not universal physics.
4. **Metacognitive self-ask:** checked whether the source is actually cold, neutrino-transparent and beta equilibrated; whether composition is present; whether the output obeys the Euler relation; whether lepton fractions close charge; whether the outer-crust correction is applicable; and whether the interpolated EOS is mechanically stable and causal between source knots.
5. **CoVe/independent derivation:** derived the first-law identities and sound-speed relation from F, derived the polynomial coefficients for PCHIP intervals and the analytic extrema conditions, and independently expanded the relativistic electron-pressure expression at small relativity parameter.
6. **Adversarial self-ask:** exercised missing composition, nonzero entropy, nonzero Q5, an equal-pressure phase boundary, acausal interpolation, muon onset, outer-crust pressure, nonstandard extra thermo columns, singular smoothing controls and a GPPVA production input rejected by the interpolation physical gate.
7. **Claim-to-evidence trace:** tied claims to source routines, test names, exact workflow runs, the manifest and the source artifact. A compile or test result is only attributed to the exact commit shown in the cited run.
8. **PDR/decision record:** retained failed CI attempts, the reason for each fix, the GPPVA refusal, unsupported capabilities, source provenance and the no-merge decision. Native subagent dispatch was unavailable and no simulated subagent review is claimed.

These labels organise reproducible work. They do not replace a derivation, source citation, test log or code diff.

## 2. Principal adversarial findings

| Finding | Severity | Response |
|---|---|---|
| A table can pass CompOSE intake and C1-C11 yet yield an inadmissible interpolated master EOS. | Critical | Added an independent master-interpolant physical gate. GPPVA_FSU2H is refused rather than silently repaired. |
| The historic global spline is not an independent physical definition. | High | Replaced it with an explicitly documented native local weighted-cubic regularisation followed by PCHIP. No exact UnivariateSpline parity is claimed. |
| A composition-free barotrope cannot determine proton, electron, muon or nuclear fractions. | Critical | Composition-free sources retain Q4 and composition placeholders; no particle fractions are fabricated. |
| Setting Q2=Q5=0 can silently change a finite-entropy or neutrino-trapped input. | Critical | Added source-row gates for negligible Q2 and Q5. |
| Electron-muon equilibrium is a scalar constraint at fixed charge fraction. | High | Added a monotone bisection solve of mu_e(Yq-Ymu)-mu_mu(Ymu)=0 and a separate muon-threshold branch. |
| A remnant nucleus with inconsistent charge closure can contaminate the core boundary. | Medium | Added an opt-in, bounded dissolution heuristic and documented that it is not a microscopic calculation. |
| The electron-pressure expression loses precision at small x_e. | Medium | Added a small-x_e series and checked its coefficients symbolically. |
| A fixed mesh is not a proof of bounds on a cubic PCHIP interval. | High | Replaced the core acceptance mesh with exact polynomial-extremum checks. The outer-crust integral bridge remains separately disclosed as mesh checked because its F is quartic. |
| A provider-specific extra column cannot be assumed to contain c_s^2. | High | Preserve by default; recompute only with explicit opt-in and Nadd=1. |

## 3. Five bounded review-and-fix cycles

### Cycle 1: First-law representation and API

**Question:** Can every returned thermodynamic field be derived from one state function?

**Change:** Added native reconstruct_master and the Python wrapper. The native routine obtains F from CompOSE Q6, optionally regularises F(log n_B), and derives all dependent thermodynamic quantities from that function. Composition-free inputs are not assigned synthetic fractions.

**Evidence:** The linear-F regression fixture tests the closed thermodynamic relations and source-array immutability.

### Cycle 2: Fortran build and outer-crust work-array correction

**Finding:** The first native build exposed a stale outer-crust electron-density reference after an array declaration had been removed.

**Change:** Replaced it with a local electron-density scalar, neval, used consistently in k_Fe and the lattice energy.

**Evidence:** The subsequent exact-head gfortran build passed.

### Cycle 3: Sample refusal, diagnostics and numerical test semantics

**Findings:** An old test still expected the function to be unported. Two analytic checks also had unnecessarily strict tolerances. The production GPPVA_FSU2H sample passes table-level C1-C11 intake but is rejected by the new master-interpolant physical gate.

**Change:** Updated the old test to require an explicit native physical rejection, not a legacy fallback; added interval and local thermodynamic diagnostic values; adjusted numerical test tolerances to cover the observed floating-point rounding residual.

**Decision:** Do not weaken physical acceptance merely to make a production table pass. The GPPVA rejection is retained as unresolved evidence that this v1 interpolant is not yet suitable for that source EOS.

### Cycle 4: Exact core-interval stability and causality

**Question:** Can a narrow unstable or acausal pocket lie between a fixed set of samples?

**Change:** Expressed each PCHIP master-function segment as a cubic in its local coordinate. Since F, mu_B and mu_B-D are cubic, while Q1 and D are quadratic, the native kernel now checks their endpoint and internal stationary-point minima. The outer-crust integral bridge is quartic in its master function, so it is checked separately on 256 points per interval and is not described as an analytic proof.

**Evidence:** The source routines cubic_polynomial_minimum and quadratic_polynomial_minimum implement the extrema; the exact-head CI exercises the full path.

### Cycle 5: Cold-source metadata and release gate

**Finding:** The algorithm sets Q2=Q5=0, so nonzero source entropy or neutrino/lepton potential must be rejected rather than overwritten.

**Change:** Added explicit Q2 and Q5 source-row gates and separate negative regression tests. Updated API, feature-parity, README and changelog documentation. Added a self-contained derivation note and this audit record.

**Evidence:** Exact-head push and PR CI both passed on be648addc3fc20a937051b5b623009c9ef2d7eb8.

## 4. Thermodynamic claim-to-evidence trace

For x=ln n_B, the implementation uses

\[
F(x)=f/n_B,\quad Q_1=F',\quad
\mu_B=F+F',\quad P=n_BF',\quad \epsilon=n_BF.
\]

It follows that

\[
D=\frac{dP}{dn_B}=F'+F'',\qquad
\frac{c_s^2}{c^2}=\frac{D}{\mu_B}.
\]

For a local PCHIP cubic F(t)=A+Bt+Ct^2+Et^3 on t in [0,1], the core gate checks:
- F(t)>0;
- Q1(t)>0;
- mu_B(t)>0;
- D(t)>0;
- mu_B(t)-D(t)>=0.

The minima of the quadratic and cubic polynomials are found from their exact derivative roots and endpoints. The outer-crust branch integrates a PCHIP Q1, making F quartic; that branch retains a dense 256-point interval gate.

| Claim | Source/test | Gate |
|---|---|---|
| One master function defines all cold thermodynamics | compose_eos_reconstruct.f90; analytic linear-F test | Formula test plus C1-C6/C10-C11 |
| Core PCHIP extrema are checked over full intervals | cubic_polynomial_minimum and quadratic_polynomial_minimum | Exact polynomial extrema |
| Composition is not invented when absent | composition-free analytic test | Exact Q4 preservation |
| Charge and baryon closure | electron-muon composition test | Yn+Yp+Yh=1 and Ye+Ymu=Yq |
| Muon equilibrium | same composition test | Independent Fermi mu_e=mu_mu check |
| Dissolution repair is opt-in and bounded | remnant-nucleus regression test | Row-level composition closure |
| Outer-crust electron/lattice pressure | outer-crust regression test | Independent pressure formula |
| Extra thermo-column semantics | extra-column regression test | Preserve by default; explicit opt-in |
| No interpolation across a first-order plateau | Iphase-boundary negative test | Explicit phase/pressure gate |
| No finite entropy or trapped-lepton state overwritten | Q2 and Q5 negative tests | Q2 approximately zero, Q5 approximately zero |

## 5. Source provenance

- S. Typel et al. (2022), CompOSE Reference Manual v3.01, [arXiv:2203.03209](https://arxiv.org/abs/2203.03209). CompOSE thermodynamic variables and table metadata.
- G. Baym, C. Pethick and P. Sutherland (1971), The Ground State of Matter at High Densities: Equation of State and Stellar Models, The Astrophysical Journal 170, 299-317, [doi:10.1086/151216](https://doi.org/10.1086/151216). Cold outer-crust nuclide sequence and Coulomb-lattice contributions.
- G. Baym, H. A. Bethe and C. J. Pethick (1971), Neutron Star Matter, Nuclear Physics A 175, 225-271, [doi:10.1016/0375-9474(71)90281-8](https://doi.org/10.1016/0375-9474(71)90281-8).
- L. Suleiman, J. L. Zdunik and P. Haensel (2024), Layers of electron captures in the crust of accreting neutron stars, Astronomy & Astrophysics 690, A301, [doi:10.1051/0004-6361/202348340](https://doi.org/10.1051/0004-6361/202348340). Modern outer-crust electron/lattice EOS treatment.

## 6. Exact-head release evidence

**Validated commit:** be648addc3fc20a937051b5b623009c9ef2d7eb8  
**Branch:** research/fortran-eos-gap-fill-v1

- Push CI: [run 42393206341](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42393206341), success.
- Pull-request CI: [run 42393211813](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42393211813), success.
- Native gfortran build: PASS.
- Frontend/direct-provider suite: 64 passed in 21.77 s.
- SHA-256 manifest: 63/63 PASS.
- Source archive generation and upload: PASS.
- Push artifact ID: 13581504022.
- Artifact size: 199,008,008 bytes.
- Uploaded artifact SHA-256: 77f40a96afd3b3f20c41742b3ca01ed1d61b6cd691162145acd271e6e95f7bf0.
- Artifact expiry: 2026-12-07.
- PR #1 remains open and unmerged; base main was not changed.

## 7. Limitations, blockers and decisions

1. GPPVA_FSU2H remains refused by this v1 interpolant's physical gate. This does not prove that the source EOS itself is nonphysical; the selected interpolation is not accepted for this input.
2. Local weighted-cubic smoothing is not legacy UnivariateSpline parity. The smoothing width, unsmoothed PCHIP and outer-crust join require source-specific convergence assessment before using an affected EOS for inference.
3. The outer-crust integral bridge has a dense mesh gate, not an analytic quartic-extremum proof.
4. The dissolution trigger is a bounded data-repair heuristic.
5. Gibbs requires two-charge surfaces P(mu_B,mu_Q), charge densities and a mixed fraction; the current 1-D input is insufficient.
6. Seidov onset stability, phase-interface TOV/Love-number jump, radial-mode stability and complete stable-branch TOV diagnostics remain unported or unverified in the native phase-transition path.
7. Microscopic RMF/Skyrme/APR/quark solvers and a complete equilibrium outer-crust nuclear-species generator remain outside this port.
8. Native subagent dispatch was unavailable. Local compilation was unavailable because the connected Desktop Commander endpoint had exhausted its monthly quota; GitHub Actions provides the exact compiler/runtime evidence.

The v1 stop condition is satisfied for the stated master-reconstruction contract. This does not close the distinct GPPVA smoothing case, Gibbs, microscopic EOS generation or phase-interface stellar-structure work.


## 8. Writing and citation-tool check

The master derivation and this audit were checked with the Academic Writing Toolkit. Both passed the British-English spelling check with zero findings. The paragraph-logic checker flagged short paragraphs around displayed equations, code examples and explicit scope exclusions; these were retained where they separate a derivation, interface contract or negative claim. The Harvard citation checker returned `notes-source-malformed` for the `**Source**` reading-note syntax emitted by its own template generator, then produced phantom-citation warnings. Its result is therefore **INCONCLUSIVE**, not PASS. Bibliographic details were checked against the CompOSE arXiv manual, the original BPS paper, the original Baym-Bethe-Pethick paper and the published A&A article. The electron-pressure low-density expansion was independently symbolically checked: the bracket is (8x_e^5/5-4x_e^7/7+x_e^9/3-5x_e^{11}/22+O(x_e^{13})), matching the native small-(x_e) series through the retained order.
