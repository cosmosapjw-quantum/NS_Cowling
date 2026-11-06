# Adversarial audit report - compose_eos v3.3.0

**Audit date:** 2026-11-06  
**Source package:** `EoS_revised_code_package.zip` (preserved; not overwritten)  
**Release scope:** original cold CompOSE EOS data, normalized data, core-crust bridge code, exact Maxwell/Gibbs paths, phase-aware TOV/Love solver, causality diagnostics, scripts, tests, documentation, and the 390-table baseline.

## Executive verdict

The v3.3 upgrade completes the previously missing *algorithmic* phase-transition paths: a sharp Maxwell construction from two locally neutral one-dimensional phase barotropes, a phenomenological CSS high-density phase constructor, a two-charge Gibbs construction from charge-resolved pressure surfaces, and a phase-aware TOV/Love solver with the pressure-corrected tidal jump. The existing 390-table full-matrix result remains **234 PASS / 156 FAIL**; the 156 FAIL cases are now attributed to the supplied APR_APR and RG_KDE0v1 high-density source tails, not to the core-crust interpolation. An APR+CSS hybrid example was generated, read back, and used in three phase-aware TOV/Love integrations. A Gibbs implementation is tested on analytic two-charge surfaces, but no Gibbs EOS is fabricated from the supplied 1-D beta-equilibrated tables because those tables do not contain the required charge-resolved surfaces.

## 1. Audit scope and evidence

The audit uses the original and normalized source trees, the 390 matched EOS outputs, the v3.3 package source, the normalizer, the complete test suite, actual saved phase-transition files, and independent recalculation of consistency and sound-speed diagnostics. Native subagent execution was unavailable in the connected environment; the work was performed as a serial audit with independent numerical checks, not presented as a parallel review.

The current release preserves raw source files under `eos_data/` and the normalized/reconstructed copies under `eos_data_normalized/`. The normalizer reproduces all 15 normalized `eos.compo` files byte-for-byte and reproduces 87/87 file-level comparisons across the normalized table files. GM_GM1's unresolved source composition remains omitted as a whole; no partial composition file is promoted.

## 2. Physical distinctions and implemented paths

### Sharp Maxwell

For two locally neutral cold phases, coexistence is defined by

$$P_H(\\mu_t)=P_Q(\\mu_t)=P_t,\\qquad \\mu_{B,H}=\\mu_{B,Q}=\\mu_t.$$

The baryon densities and energy densities may jump. At fixed `P_t` and `mu_t`, thermodynamic consistency gives

$$\\Delta\\epsilon=\\mu_t\\Delta n_B.$$

The new `build_maxwell_eos` implementation solves the pressure crossing in `P(mu_B)`, requires monotone stable causal input phases, and inserts adjacent rows at common pressure and chemical potential with `n_Q > n_H` and `epsilon_Q > epsilon_H`. Readback verifies the table values against the transition sidecar. This is distinct from the legacy `method="maxwell"` core-crust bridge, which remains a smooth compatibility interpolant and is **not** an exact Maxwell transition.

The optional CSS constructor defines a high-density branch with a supplied `Delta epsilon` and `c_s^2`:

$$\\mu(n)=\\mu_t(n/n_Q)^s,\\quad P(n)=P_t+\\frac{s}{1+s}(\\epsilon_Q+P_t)[(n/n_Q)^{1+s}-1],\\quad \\epsilon=\\mu n-P.$$

This is a phenomenological constant-sound-speed phase. It is not a substitute for a microscopic quark-matter calculation. The high-phase charge/composition profile is not reconstructed; the output marks composition unavailable and records the effective one-fluid assumption.

### Gibbs mixed phase

The Gibbs path requires both phases as two-charge pressure surfaces `P(mu_B,mu_Q)`, with derivatives `n_B=partial P/partial mu_B` and `q=partial P/partial mu_Q`, plus a lepton sector. It solves pressure equality and global charge neutrality:

$$P_H(\\mu_B,\\mu_Q)=P_Q(\\mu_B,\\mu_Q),$$

$$ (1-\\chi)q_H+\\chi q_Q+q_l=0. $$

The mixed baryon density is `(1-chi)n_H + chi n_Q`. The supplied package tables are one-dimensional beta-equilibrium barotropes. They do not contain `P(mu_B,mu_Q)` surfaces or their charge susceptibilities. The Gibbs API is therefore tested against analytic synthetic surfaces, but deliberately does not claim to produce a Gibbs EOS from the supplied one-dimensional data.

### TOV and tidal jump

The phase-aware solver integrates the high-density phase to `P_t`, applies the `l=2` logarithmic-derivative jump, then continues through the low-density phase. The implementation uses the corrected finite-pressure condition of Takátsy and Kovács (2020):

$$y(r_t^+)-y(r_t^-)=-\\frac{4\\pi r_t^3\\Delta\\epsilon}{m(r_t)+4\\pi r_t^3P_t},\\qquad G=c=1.$$

The interface record stores `r_t`, `m_t`, geometric `P_t`, geometric `Delta epsilon`, the two `y` values, and `Delta y`. The saved TOV result's `Delta y` was independently recomputed from those stored interface quantities; the residual is at floating-point roundoff. The solver reports the Seidov small-core onset condition separately and does not claim a complete radial-mode stability analysis.

## 3. Original 390-table FAIL attribution

The independent failure classifier parsed all 390 matched EOS outputs and reran the hard gates. Results:

| Status/root cause | Count |
|---|---:|
| PASS | 234 |
| SOURCE_CORE_ACAUSALITY_PRESERVED | 156 |
| Matching/numerical candidate | 0 |
| Audit read errors | 0 |

The 156 failures split into 78 APR_APR and 78 RG_KDE0v1 outputs. Source profiles calculated with PCHIP, natural cubic spline, and centered finite differences all show `c_s^2 > 1` in the same high-density tails. For all 156 failing outputs, the post-transition matched arrays `(n_B,P,epsilon)` are unchanged from the corresponding source core tail within the recorded tolerance. The full derivation and numerical table are in `FAILURE_ROOT_CAUSE_REPORT_v3_3.md`.

This is a source-model limitation over the full tabulated high-density domain, not evidence that the true physical EOS must be acausal. It is also not a cross-EOS matching conflict. The different crust/core saturation properties remain a possible physical systematic in the junction location and structure, but they do not explain these C10 failures.

## 4. Causal-completion experiments

The optional completion replaces the superluminal source tail with a thermodynamically reconstructed constant-`c_s^2` continuation, preserving `P`, `epsilon`, and `mu_B` at the chosen join and leaving the original input untouched.

- APR_APR, tail `c_s^2=0.95`: whole consistency gate PASS; three independent derivative maxima below unity (PCHIP 0.9818, natural cubic 0.9953, centered FD 0.9809). This is a robust alternative EOS under the tested diagnostics, not the original APR table.
- RG_KDE0v1, tail `c_s^2=0.80`: package hard gate PASS, but independent PCHIP/cubic-spline diagnostics remain above unity near the join. It is retained as diagnostic-only; it is not a robust promoted alternative.

The code's CLI now emits both full consistency and cross-scheme diagnostics and labels an output `diagnostic_only` unless all checks pass.

## 5. Actual APR+CSS Maxwell/TOV validation

A phenomenological hybrid was constructed from the APR_APR low-density branch at `n_t=0.60 fm^-3`, `P_t=123.020157 MeV fm^-3`, `Delta epsilon=120 MeV fm^-3`, and `c_s,Q^2=0.75`. The APR hadronic branch sound speed at the chosen point is approximately 0.566. The read-back Maxwell table passes the full consistency check.

The transition is:

- `mu_t = 1242.932120 MeV`;
- `n_H = 0.600000 fm^-3`, `n_Q = 0.696546 fm^-3`;
- `epsilon_H = 622.739115 MeV fm^-3`, `epsilon_Q = 742.739115 MeV fm^-3`;
- `Delta epsilon / Delta epsilon_critical = 0.241984`, below the Seidov small-core threshold for this parameter set.

For a central pressure of `815.642916 MeV fm^-3`, the tightest of three integrations gives approximately `M=1.96925461 Msun`, `R=9.762497 km`, `k2=0.01838394`, `Lambda=5.227468`. These are the outputs of the stated phenomenological CSS model, not predictions for the original APR EOS. Across `(rtol,max_step)=(5e-7,0.15 km)`, `(1e-8,0.08 km)`, and `(5e-9,0.04 km)`, the loose-to-tight relative differences were about `3.6e-8` in mass, `7.6e-4` in radius, `3.8e-3` in k2, and `5.7e-5` in Lambda. The TOV convergence JSON preserves all runs. The shipped CLI path `construct_css_maxwell_eos.py` -> `run_phase_transition_tov.py` was independently executed and returned SUCCESS; `causal_complete_eos.py` also wrote a post-readback report with `diagnostic_only=false` for the APR 0.95 alternative. Both installed wheel console entry points (`compose-eos`, `compose-tov-audit`) returned help successfully.

## 6. Five review/fix rounds

| Round | Scope | Adversarial check / change | Outcome |
|---|---|---|---|
| R1 | Source and matrix discovery | Recomputed source causality by three derivative schemes and compared all 156 matched high-density tails to their core inputs | Established source-tail cause; no matching-candidate failures |
| R2 | Step-back and physical API design | Separated 1-D locally neutral Maxwell phases from 2-D charge-resolved Gibbs phases; implemented exact Maxwell, CSS, Gibbs, phase-aware TOV and jump condition | Algorithms implemented with explicit input contracts |
| R3 | Test-driven adversarial debugging | Initial TOV test exposed geometric-vs-CompOSE pressure-unit misuse; Maxwell test exposed high-energy spline/scalar name collision | Both defects repaired; focused tests pass |
| R4 | Thermodynamic and provenance hardening | Required phase source gates and Euler identity; preserved charge-potential fields; validated sidecars; added Gibbs convexity check and Gibbs roundtrip; challenged KDE causal-tail with three schemes | Focused phase tests pass; KDE causal alternative remains diagnostic-only |
| R5 | End-to-end release verification | Actual APR+CSS Maxwell construction, sidecar readback, independent `Delta y` reconstruction, three TOV convergence runs, full matrix failure classification, full test suite, compile/wheel/ZIP checks | Release artifacts generated with explicit limitations and hashes |

The named audit sequence was executed as an evidence workflow, not as a claim to expose internal reasoning: self-discover (inspect actual source/outputs), step-back (classify available thermodynamic state variables), CRAG (retrieve and corroborate primary sources), metacognitive self-ask (list assumptions and unavailable inputs), CoVe (independent derivative and symbolic checks), adversarial self-ask (challenge unit, phase-order, jump-sign and provenance assumptions), CCoT (claim-to-equation-to-code-to-test lineage), and PDR (plan/do/review with five bounded correction rounds). Native subagent dispatch was not available and was not simulated.

## 7. Verification summary

- Full regression suite: 54 passed (one non-fatal Matplotlib `Axes3D` import warning unrelated to EOS calculations).
- Phase-transition focused suite: 7 passed.
- Python `compileall`: passed.
- Wheel build and isolated import/CLI smoke test: recorded in the delivery manifest.
- 390 matched tables: independently read and gated; 234 PASS / 156 source-core C10 FAIL.
- APR CSS Maxwell: CompOSE write/read, transition sidecar checks and full consistency gate PASS.
- Phase-aware TOV: three tolerance/max-step runs; independent tidal jump formula residual at roundoff.
- Gibbs: analytic 2-charge surface construction, global charge neutrality, pressure equality, convexity checks and write/read sidecar test.

## 8. Remaining limitations

1. The package's original 156 APR/KDE0v1 full-table FAILs remain FAIL. Only optional alternative causal completions are provided.
2. The KDE0v1 completion at `c_s^2=0.80` passes the package hard gate but remains method-sensitive; it is not promoted as robust causal.
3. Gibbs construction cannot be run on the supplied 1-D beta-equilibrated EOS tables; actual charge-resolved two-dimensional phase surfaces are required.
4. The CSS hybrid high-density phase is phenomenological. Its charge/composition profile is not microscopically reconstructed.
5. The TOV example is one fixed-central-pressure solution, not a complete mass-radius sequence or 390-model TOV campaign. Seidov is a small-core onset criterion, not a full radial stability computation.
6. The TOV surface is the first pressure point in the tabulated EOS. No arbitrary extrapolation to vacuum is claimed.

## 9. References

- CompOSE Reference Manual v3.00, Sec. 4.2.3, Eqs. (4.20)-(4.22): https://compose.obspm.fr/download/pdf/manual_v3.00.pdf
- Takátsy & Kovács, Phys. Rev. D 102, 028501 (2020), corrected finite-pressure tidal jump: https://doi.org/10.1103/PhysRevD.102.028501; OA text: http://real.mtak.hu/114802/1/Postnikov_comment_v2.pdf
- Postnikov, Prakash & Lattimer, Phys. Rev. D 82, 024016 (2010): https://doi.org/10.1103/PhysRevD.82.024016
- Bhattacharyya, Mishustin & Greiner, Maxwell vs Gibbs construction, J. Phys. G 37, 025201 (2010): https://arxiv.org/abs/0905.0352
- Han & Steiner, tidal deformability with sharp phase transitions, Phys. Rev. D 99, 083014 (2019): https://arxiv.org/abs/1810.10967
- Fortin et al., unified crust-core EOS and systematic effects, Phys. Rev. C 94, 035804 (2016): https://doi.org/10.1103/PhysRevC.94.035804
- Lindblom, causal representations of neutron-star EOS, Phys. Rev. D 97, 123019 (2018): https://scispace.com/papers/causal-representations-of-neutron-star-equations-of-state-27ptzdab5t
