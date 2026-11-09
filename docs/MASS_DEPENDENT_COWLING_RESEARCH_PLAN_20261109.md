# Mass-dependent crust-core matching response: research plan

Date: 2026-11-09  
Track: historical smooth crust-core matching cohort; physically separate from native sharp-Maxwell branch  
Status: implementation staged on a separate research branch; exact-head pilot CI is the next acceptance gate

## 1. Scientific objective

Determine how thermodynamically admitted crust-core matching prescriptions alter neutron-star structure and the relativistic Cowling f- and p1-mode spectrum as a function of stellar mass along the initial stable TOV branch.

The central question is:

> How does thermodynamically admissible crust-core matching uncertainty propagate through neutron-star structure and non-radial oscillation spectra across the stellar mass sequence, and how does that response depend on core and crust families?

The intended contribution is not another single-mass comparison. It is a controlled, paired, crossed-cohort map of mass-dependent response and its core/crust dependence, with explicit numerical and mode-identity gates.

## 2. Frozen input authority and scope

The new execution branch is research/mass-dependent-cowling-v1, forked from live native-Maxwell feature head a09ec30d74b363652d2d1aefe990557c468679dd. PR #3 remains open, draft and unmerged; this track does not modify it.

The frozen smooth-cohort authority is:
- oscillation_release_v1/data/spectrum_summary_234.csv (234 rows; target mass 1.3300106263 solar masses);
- oscillation_release_v1/data/eos_projected/ (projected tables);
- oscillation_release_v1/provenance/release_summary.json;
- docs/V5_PAIRED_MATCHING_SENSITIVITY_20261109.md;
- docs/EXISTING_SPECTRUM_REUSE_VERIFICATION.json.

The crossed cohort is 13 crusts × 3 admitted core families × 6 smooth matching prescriptions: 39 fixed core/crust pairs and 234 method rows. The three cores are BBB_BHF2, GM_GM1 and RG_SLY4. The six methods are css, cubic_F, cubic_P, maxwell, quintic_F and quintic_P. These are not 234 independent microscopic EOS draws.

The historical method token maxwell is a smooth near-pressure-plateau bridge. It has no density discontinuity or latent heat and is not a first-order Maxwell phase transition. The native sharp-Maxwell solver is a separate physical track with distinct interface transfer conditions and conversion assumptions.

The prior spectrum-reuse audit reports 234 rows and recomputed_modes=0. All old mode frequencies are fixed at approximately 1.3300106 solar masses. They are a frozen reference anchor, not mass-dependent evidence.

## 3. Hypotheses and falsifiers

H1: low-mass p1 amplification. The paired six-method p1 range may increase toward low masses and contract at higher masses. Canullán-Pascual et al. (2025) already report low-mass sensitivity, so reproducing that qualitative behaviour is an external-consistency check, not novelty by itself.

H2: core-dependent response shape. The mass where the six-method range peaks, the width of the high-sensitivity interval, or the sign pattern of method deviations may depend on the core family. This is a candidate differentiator, not an assumed result. It is falsified if core-conditioned curves collapse within numerical uncertainty.

H3: crust-conditioned response. At fixed core, the method response may depend on crust family and may be concentrated in particular core/crust pairs rather than follow a universal mass-only curve. It is falsified if the response is consistent with a common mass-only curve.

H4: mode contrast. The p1-mode may remain more matching-sensitive than the f-mode across a substantial part of the stable branch. It is falsified if the contrast disappears after mode-identity and convergence gates or is confined to a small subset.

H5: physical mechanism. A later analysis may connect response to crustal acoustic travel time and mode amplitude near the matching region. No sensitivity-kernel result is assumed. A kernel claim requires a variational derivation for the relativistic Cowling system, including boundary and matching terms.

## 4. Confirmatory mass design and stable branch

The predeclared target masses are
\[
M/M_\odot = 0.8,\ 1.0,\ 1.2,\ 1.33,\ 1.4,\ 1.6,\ 1.8,\ 2.0.
\]

For every EOS:
1. Scan a logarithmic central-pressure grid below the projected EOS pressure ceiling, anchored to the frozen 1.33-solar-mass central pressure.
2. Integrate TOV at every sample. Preserve failures explicitly.
3. Extract only the initial connected rising-M(Pc) branch. The current implementation resolves a first maximum after two consecutive fractional mass decreases larger than 1e-4. A failed sample is a hard barrier; later points are not joined across it.
4. If the scan terminates before a maximum is resolved, report a branch endpoint, not Mmax.
5. For each mass bracketed inside the branch, pass the adjacent sampled central pressures to the production Fortran target-mass solve.
6. Do not extrapolate a target beyond the resolved branch or across a failed-sample barrier.

The turning-point criterion is a barotropic branch/stability proxy. Unqualified radial stability claims require a separate radial-pulsation calculation, especially for genuine phase transitions or disconnected branches.

## 5. Spectrum solve and mode identity

The frozen f and p1 frequencies seed search windows using only the mean-density scaling
\[
\nu_{\rm dyn}\propto\sqrt{M/R^3}.
\]
The predeclared window is 0.72–1.28 times the scaled reference frequency. This is a root-bracketing device, not a fitted physical prediction.

The `cell_solve` EOS-cell integrator now accepts a configurable `r_limit_m` in its namelist; its historical default remains 1e6 m. For this mass scan, the explicit exploratory ceiling is 1e8 m and is recorded in every model manifest. The wrapper adds `solve_cowling_tov_star(..., task=star)` so sequence scans use the same pressure-cell integration and central-start contract as the mode calculation. The generic `case_driver real_table_star` API is retained separately. The Python wrapper accepts optional explicit pc_bracket_pa=(Pc_lo,Pc_hi). When omitted, its former 0.65–1.45 multiplicative bracket remains unchanged. When supplied, the bracket is validated and passed into native cell_solve. This constrains the target-mass search to the sampled rising branch.

A spectrum row is admitted only if:
- the target-mass residual closes;
- both shooting residuals satisfy the declared tolerance;
- the smooth l=2 Cowling W/V node signature is established in both profile segments:
  - f: zero W and V internal nodes in inner and outer segments;
  - p1: one W and V internal node in the inner segment, zero in the outer segment.

The identity contract is smooth_l2_WV_nodes_v1. It is not transferable to the sharp-Maxwell interface branch.

## 6. Paired estimands

For y in {f,p1}, at fixed core c, crust s, mass M and method m, define
\[
\Delta_m^{(y)}(M;c,s)
=100\left[\frac{y_m(M;c,s)}{\bar y_6(M;c,s)}-1\right],
\qquad
\bar y_6=\frac{1}{6}\sum_{m=1}^{6}y_m.
\]

The six-method paired range is
\[
{\cal R}_6^{(y)}(M;c,s)
=100\,\frac{\max_m y_m-\min_m y_m}{\bar y_6}.
\]

An estimand enters a paired summary only if all six methods for that fixed core/crust/mass cell have valid, identified roots. Missing methods are not imputed. At every mass report the complete-pair denominator, median, 90th percentile and maximum range, plus core/crust-resolved maps and method deviations.

The 39 fixed core/crust pairs are the physical pairing units. Do not treat the 234 method rows as independent EOS samples, posterior samples or population-level draws. The existing 1.33-solar-mass variance decomposition remains descriptive for its selected finite cohort, not causal or posterior.

## 7. Numerical validation gates

N1: frozen baseline reproduction. Recompute TOV M,R at every reference central pressure using the same `cell_solve task=star` EOS-cell integrator as the spectrum calculation, with the recorded v5 settings `r0_m=10 m`, `rtol=1e-9`, `step_cap=0.05` and the identical first-row finite surface-pressure cutoff. The separate `case_driver real_table_star` path is not an acceptable direct reproduction comparator because its ODE and central-start contracts differ. Current pilot acceptance requires relative M,R differences no larger than 2e-5.

N2: target-mass solve. Require relative mass residual no larger than five times the configured mass_rtol; default mass_rtol is 1e-10.

N3: shooting and identity. Require both matching residuals no larger than 1e-6, plus the exact W/V node signature above. A small determinant alone is insufficient.

N4: branch resolution. Preserve all pressure samples and failures; state whether the first turning point was resolved. A truncated endpoint is never relabelled Mmax.

N5: refinement versus effect size. For representative low-mass tail, 1.33-solar-mass anchor and high-mass cases, vary TOV/Cowling tolerance, maximum step and root brackets independently. Frequency movement must be at most one tenth of the smallest paired method effect reported. Otherwise weaken or withhold the corresponding physical claim.

N6: independent benchmark. Compare representative EOS/mass points with published relativistic Cowling/full-GR results or a separately maintained solver. Cowling/full-GR systematic offsets are distinct from within-Cowling matching shifts and must be tabulated separately.

N7: recognized matching prescriptions. The six historical smooth bridges are not a sufficient substitute for direct literature baselines. Add (i) simple EOS concatenation as a deliberately non-thermodynamic diagnostic control only, (ii) linear P(epsilon) matching, (iii) the thermodynamically consistent P(epsilon) interpolation used by Canullán-Pascual et al. (2025), reconstructed from the published prescription, and (iv) a CUTER-style or otherwise documented thermodynamically consistent crust reconstruction. For every method preserve endpoint data, units, source checksums, chemical-potential/first-law closure, causality and interpolation receipts. Do not label simple concatenation physically admissible, and do not call an approximate CUTER-like implementation CUTER itself.

## 8. Structural and mechanism extensions

This increment produces TOV M(Pc), branch endpoints/turning-point candidates, R(M), and f(M),p1(M). It does not yet produce crust mass, crust thickness, tidal Love number/Lambda, full-GR complex mode frequencies or damping times.

Next mechanism work:
- extract the transition radius and enclosed mass from pressure-event output or validated profile interpolation;
- compute crust acoustic travel time with the correct relativistic Cowling radial measure and EOS sound speed;
- compare normalized eigenfunctions and overlap near the matching region;
- derive a variational/sensitivity kernel only after all boundary and matching terms are included and checked.

Next astrophysical work:
- implement Love-number/tidal deformability and source-backed observational likelihoods;
- combine NICER mass-radius posteriors, maximum-mass and tidal constraints;
- admit CompOSE/source-backed EOSs only with source checksums, units, phase labels and thermodynamic admission receipts.

## 9. Publication decision

JKPS / EPJA: viable only after the full crossed cohort has mass-dependent paired spectra, explicit stable-branch coverage, refinement below one tenth of the smallest reported effect, and independent benchmarks. The distinguishing result must be a reproducible core/crust-conditioned mass-response structure, not simply “p1 is more sensitive than f”.

A&A / MNRAS: requires realistic EOS provenance, structural observables, tidal deformability and observational inference.

PRD: requires a new relativistic perturbation/sensitivity derivation or systematic Cowling/full-GR response result.

## 10. Implementation status

Stable-branch extraction, adjacent-Pc bracketing, six-method complete-pair estimands, explicit Pc and radius-limit controls, the shared cell_solve TOV-only backend, the all-cohort runner and the pilot/full workflow are staged on the new branch. The first integration run (Actions 42598463512) exposed a non-comparable baseline: the runner used `case_driver real_table_star` (`r0=1 m`, `rtol=1e-8`) while the v5 reference metadata records `cell_solve` (`r0=10 m`, `rtol=1e-9`, `step_cap=0.05`). The same run also reached the old hard-coded 1e6-m direct-cell radius ceiling. The wrapper now uses the same `cell_solve task=star` backend and exposes the ceiling without changing its default. The validator syntax defect is corrected. A second run (42599325312) passed the four helper tests but failed at a Python keyword/signature mismatch before native compile or TOV integration; the corrected `solve_model` signature is the next patch. No mass-dependent spectra or Mmax result has passed yet.

References:
1. Canullán-Pascual et al. (2025), Astronomische Nachrichten, DOI 10.1002/asna.20240150, arXiv:2502.02373.
2. F. Nola (2026), A&A 713, A34, DOI 10.1051/0004-6361/202659810.
3. Kunjipurayil et al. (2022), Phys. Rev. D 106, 063005, arXiv:2205.02081.
4. Pradhan and Chatterjee (2022), Phys. Rev. C 106, 015805, DOI 10.1103/PhysRevC.106.015805.


## Provenance qualification for the v5 backend

The current solve_v5_* adapter compiles the 35 source paths listed, in their recorded order, by oscillation_release_v1/workflow/compile_manifest.json and records a runtime build key over those branch source bytes plus the manifest. The release_summary.json stores a reference source_zip_sha256 of b2c0dcaf6a5b07956cf8d33c39651d4ade32f44b4c49b33aa08fd84fb458185b. The Dropbox source archive was not downloaded and independently hashed in this session, so the adapter is manifest-path/order matched, not yet byte-identity certified against the archived release source. The next CI artifact must retain the runtime build key and manifest SHA.


## 11. Forward-sensitivity route for the mechanism analysis

The current smooth Cowling W,V system has been recast as y'=A y. The parameter tangent S_zeta'=A S_zeta+A_zeta y, together with the frequency tangent, yields the simple-root derivative -D_zeta/D_omega. Composition with the fixed-mass identity gives a computable response at constant gravitational mass. The derivation and a synthetic ODE finite-difference check are recorded in docs/COWLING_FIRST_VARIATION_RESPONSE_20261109.md.

This is an implementation route, not yet a physical sensitivity kernel. A physical calculation must differentiate the TOV background, EOS-cell closure, centre expansion, finite-positive-pressure surface/event, matching determinant and mode identity. For a sharp Maxwell interface, the slow/rapid transfer map must be differentiated separately. The existing sharp-interface track is not silently folded into the smooth six-method cohort.

The q-path mock was strengthened from mesh-only admission to whole-polynomial-interval extremum checks. Its exact lower bound applies only to the synthetic path and does not establish that the six historical matching methods lie on one admissible differentiable path.


## 12. Branch topology, radial stability and the secondary-branch stop rule

The first connected rising-M(Pc) segment remains the conservative mass-root search proxy, but it is not itself a radial eigenvalue calculation. The first 24-point exact-v5 pilot shows a first mass maximum near 0.7432 Msun and R near 1232 km, followed by a large excursion and a later rising interval reaching about 1.915 Msun at the scan endpoint. This same topology occurs across all six methods in the pilot pair, but the coarse grid does not establish whether it is a physical transition structure, an unresolved numerical feature, or an unstable/stable branch sequence.

The code increment adds a kinematic monotonic-segment partition and raises the six-model pilot grid from 24 to 96 Pc samples. It labels all candidate segments as radially unresolved. No secondary branch is sent to the Cowling mode solver.

A direct radial-mode calculation is required before secondary branches are admitted. In geometric units, the standard radial problem uses the dimensionless Lagrangian displacement xi and DeltaP:
- dxi/dr = (nu' - 3/r)xi - DeltaP/(r Gamma_1 P);
- dDeltaP/dr = {exp(2lambda)[omega^2 exp(-2nu)-8 pi P] + nu'(4/r+nu')} (epsilon+P)r xi - {nu'+4 pi(epsilon+P)r exp(2lambda)}DeltaP.
At the centre, DeltaP(0)=-3 Gamma_1 P_c xi(0); at a zero-pressure surface DeltaP(R)=0. The current v5 background stops at finite positive pressure, so an outer-EOS continuation or a separately derived finite-cutoff boundary condition is required.

The closure is part of the physics. Gamma_eq=(epsilon+P)/P*dP/d epsilon tracks the barotropic equilibrium EOS; Gamma_1=(epsilon+P)/P*(partial P/partial epsilon)_(s,Y_i) describes an adiabatic perturbation at frozen entropy/composition. A radial-stability result must state which closure was used. For a genuine sharp interface, the perturbative junction condition is independent of the equilibrium interpolation: slow and rapid conversion have different [xi,DeltaP] transfer laws. The historical six methods remain smooth bridges and receive no sharp-interface jump condition.

Full equations and the exact 24-point pilot evidence are recorded in docs/MASS_SEQUENCE_BRANCH_STABILITY_AUDIT_20261109.md. The next exit sequence is: 96-point TOV topology -> 96-to-192 refinement -> same-interpolation independent TOV reference -> direct radial omega_0^2 solver and closure audit -> only then target-mass f/p1 spectra on a certified branch.


## 13. 96-to-192 topology refinement

The first sampled turning point moved from M about 0.7433 Msun, Pc about 0.172 MeV fm^-3, R about 1233 km at 24 samples to M about 0.7651 Msun, Pc about 0.227 MeV fm^-3, R about 1880 km at 96 samples. This is not converged. A 192-point run is now the pilot configuration. The maxwell compatibility bridge has a distinct extra turning pair inside its low-sound-speed interval; all such pieces remain radial-stability unresolved.


## 14. 192-to-384 refinement and independent TOV cross-check

At 192 Pc samples, all six rows completed with no failed TOV integrations and v5 M/R anchors passed. The first maximum mass moved by only about 2.2e-4 relative from the 96-point estimate, but its radius moved by about 1.6%, while the first subsequent minimum radius shifted from about 418 km to about 195 km in the css row. This is insufficient convergence. The workflow is moved to 384 Pc samples.

The next CI increment also includes an independently coded SciPy DOP853 direct-r TOV integrator. It reconstructs the exact piecewise cubic-Hermite P(mu) table from raw native knots, using the raw baryon density n as the endpoint derivative dP/dmu. It checks five central pressures spanning the first fold, falling segment, the extra maxwell smooth-bridge folds, and a high-density reference against frozen v5 cell_solve task=star. The target is <=1e-6 relative M and R discrepancy. This verifies structure integration and EOS interpolation equivalence only. It does not classify radial stability or validate f/p1 modes.


## 15. 384-point sequence and independent structural validation

The six-model 384-point native TOV scan completed with 2304/2304 valid integrations and all frozen-v5 M/R anchor checks passing. The first sampled maximum mass is now 0.7652866-0.7652926 Msun, stable at about 2.6e-5 relative compared with 192 points. Its radius still moves by about 0.75%, and the first post-maximum minimum radius remains non-monotonic across 96/192/384. The low-density turning topology is therefore only partially converged.

A direct-r SciPy DOP853 calculation using the same raw-knot cubic-Hermite P(mu) EOS was compared with the Fortran q-coordinate integrator. The initial coarse comparison gave maximum relative mass difference 3.60164e-7 and radius difference 1.46572e-6; the predeclared 1e-6 radius criterion was narrowly missed. The acceptance was not relaxed. The new repair stores casewise output on failure, tightens both integrators and retains the same 1e-6 criterion for the refined comparison.

Next: inspect the per-case refined TOV receipt; investigate any remaining mismatch through the EOS state and surface-event residuals. Only after the structural comparison is accepted should the selected M(Pc) folds be fitted/refined locally. Direct radial omega_0^2, an explicit Gamma_1 closure and a physically valid surface condition still gate all f/p1 mass-dependent spectra.


## 16. Exact-head independent TOV acceptance and physical next gate

The exact-head run 42617251360 completed the 384-point six-model sequence (2304/2304 valid TOV integrations) and all six frozen-v5 M/R anchor reproductions. The same-EOS independent SciPy DOP853 check produced a persisted PASS receipt at five central pressures. Its maximum tight SciPy/native relative errors were 1.6233e-8 in mass and 1.1433e-7 in radius, below the predeclared 1e-6 criterion; the maximum native refinement changes were 1.30e-10 in mass and 7.06e-11 in radius. This closes only the selected-point structural TOV comparison, using identical raw-knot cubic-Hermite P(mu) interpolation and the shared finite-positive-pressure surface cutoff.

The Actions run still reports failure because the final CLI summary used obsolete v1 key names after the v2 receipt had been written. The code now maps the summary directly to the v2 receipt schema, with a regression test. This does not relax the numerical criterion or the physical validator.

The physical sequence gate remains blocked: 0/48 f/p1 roots were admitted, no complete six-method target-mass cell exists, and the first rising-proxy endpoint is about 0.7653 Msun. The later high-density rising segment is only a kinematic candidate; it is not declared radially stable. The next substantive physics task is a direct radial pulsation calculation for the nodeless fundamental eigenvalue omega_0^2, after fixing the closure (barotropic-equilibrium versus frozen-composition Gamma_1) and deriving/implementing a valid boundary condition for the finite-positive-pressure surface. If composition derivatives are absent, restrict claims explicitly to barotropic-equilibrium stability. Only then attempt f/p1 roots on branches whose radial eigenvalue and node count are validated.


## 17. Finite-pressure surface continuation before radial eigenmodes (historical R13 diagnostic; superseded for primary analysis)

**Supersession:** the synthetic C1 tangent/curvature tail described in this historical section is retained as a mathematical diagnostic only. It is not selected as a physical crust and is not used in the current primary analysis. R14 replaces this work path with exact native pressure-knot truncation of the unmodified projected EOS tables. No new EOS rows are created.

The projected EOS is terminated at (P_{cut}=9.000000079422\times10^{-12}\;\mathrm{MeV\,fm^{-3}}), with (\epsilon_{cut}=4.354134734530827\times10^{-7}\;\mathrm{MeV\,fm^{-3}}), (n_{b,cut}=4.67962\times10^{-10}\;\mathrm{fm^{-3}}), and (\mu_{b,cut}=930.4654511545\;\mathrm{MeV}). The first five raw knots are exactly identical across the six matching-method EOS files, so a shared outer continuation is available without introducing a method-specific exterior difference.

The right-sided Hermite estimate gives (c_{s,cut}^2=3.1507096132423\times10^{-5}) and (\Gamma_{eq,cut}=1.52432195647). The local right-cell (\Gamma_{eq}) decreases to approximately 1.46155 by the fifth knot. Log-log fits of (P) versus (\epsilon) over 5, 10, 20, 40, 80 and 160 knots give exponents 1.58002, 1.54586, 1.46579, 1.39937, 1.36705 and 1.32284, respectively, with increasing log-residual. A single broad-window polytrope is therefore not a robust extrapolation rule.

Historical R13 reference model only (not selected for the primary analysis):

The former diagnostic used a C1 tangent relation P=K epsilon^gamma_s with gamma_s=(epsilon_cut/P_cut)c_s,cut^2 and completed n_b,mu_b through the cold first law. Its first-law algebra is a valid barotropic closure, but it supplies no microphysical outer-crust composition and is an unsourced extrapolation below the input table. The corresponding synthetic-tail workflow is manual-dispatch only and is not part of the EOS-preserving primary path.

R14 instead terminates the unchanged source EOS at exact native pressure knots and reports finite-pressure internal-isobar M(r_cut), r_cut. This measures the truncation response only within the provided positive-pressure domain; it cannot determine the missing physical crust or P=0 radius.

The first radial pulsation calculation may set (\Gamma_1=\Gamma_{eq}), but results must be called barotropic-equilibrium radial stability. Frozen-composition (\Gamma_1) is not available from the current projected table without composition derivatives. Only after a sourced physical low-density EOS and a justified P=0 surface boundary have been established should the code solve for the nodeless fundamental (\omega_0^2), check node count and compare against published radial-mode reference models. Mass-dependent f/p1 calculations remain blocked until that radial eigenvalue gate is complete.


## 18. Exact-head verification after the CLI repair

Actions run 42618432567 (source-code HEAD `2e4dc9982238bd937a92019dd46580b37923998f`) completed all computational steps. It passed 10 mass-sequence tests, 6 fixed-mass/Cowling mock tests and 5 Hermite-EOS/CLI tests. The six-model TOV scan returned 2304/2304 valid points, no model failures, and no frozen-v5 M/R baseline failures. The independent same-EOS SciPy DOP853 comparison passed its unchanged (10^{-6}) criterion; maximum tight relative mass/radius discrepancies were (2.1897\times10^{-7}) and (1.3524\times10^{-7}). The CLI receipt-schema KeyError is fixed and covered by regression test.

The workflow correctly stops at the physical validator: the first rising-branch proxy reaches only about 0.7653 Msun, below the required 1.0 Msun mass; all 48 requested spectrum rows are outside that proxy, so no complete six-method paired cell exists and no f/p1 spectrum is admitted. Overall workflow conclusion is failure by design at the unresolved science gate, not a numerical-integrator or CLI-code failure. Artifact ID 13675450289 and SHA-256 `30922c2fafd5260d3c03d6852057530555cd997a64962d5b59d3f8af966af2ee` were verified after download.

R14 has replaced the synthetic continuation as the active primary surface-boundary investigation. The active result is an EOS-preserving scan at eight existing native pressure knots for six matching methods (48 TOV integrations), with exact source-table hashes checked before and after and native background profiles exported. All 48 cases passed; source tables remained unchanged; the six native-edge M/R baselines reproduced within 2e-5. This establishes truncation convergence within the tabulated positive-pressure domain only. It does not establish a P=0 surface or a complete-star radius. The next physical gate is a sourced low-density crust EOS, or a separately specified finite-domain boundary-value problem whose truncated-domain interpretation is explicit. Radial \\(\\omega_0^2\\) and full-star Cowling f/p1 calculations remain blocked until a justified surface condition is available. Initial (Gamma_1=Gamma_{eq}) would support only a barotropic-equilibrium stability claim, not frozen-composition stability.

## 19. R14 native-table-only boundary scan (2026-11-09)

### Authority and invariants

The six projected EOS CSVs in the frozen GPPVA_FSU2 x BBB_BHF2 cohort remain the sole primary EOS inputs. Each contains 323 native rows. Their Git blob identities and before/after SHA-256 values are recorded in the R14 artifact. No row was appended, deleted, rewritten, sorted, smoothed, clipped, refitted, or replaced. The frozen v5 native-knot cubic-Hermite P(mu_b) provider is used as the continuous representation between supplied knots; this does not create or write a new EOS table.

### Design

The selected pressure cutoffs are nearest existing rows 1, 3, 4, 6, 7, 11, 15 and 20, corresponding to actual log10 offsets 0, 0.311229, 0.467474, 0.779805, 0.936080, 1.560687, 1.994697 and 3.049259 dex above the first native row. Each method is run at the same frozen reference central pressure with the existing v5 TOV backend and settings r0=10 m, r_limit=1e8 m, rtol=1e-9 and step_cap=0.05. The native background CSV is exported at every finite-pressure endpoint.

### Exact-head verification

Actions [42640660436](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42640660436) passed on exact HEAD `aff4e6ad2789555648f2e564a96bbabd952ee4ac`.

- 48/48 TOV cases completed; zero failures.
- 6/6 original EOS source CSVs retained their exact SHA-256 before/after values and expected Git blob identities.
- 6/6 row-1 table-edge baseline reproductions passed the unchanged 2e-5 relative M/R tolerance.
- 48/48 backend-exported background profiles passed SHA-256 and row-count validation.
- Zero synthetic EOS tables were emitted.
- Artifact ID `13684872040`, 18,300,545 bytes, SHA-256 `e124339e932636e9e08c92298b2e1269aca3fac4980230083a28a3451421d8f3`.

At source row 20 relative to row 1, the enclosed mass inside the selected isobar changes by approximately (-4.54\times10^{-10}), while the isobar radius decreases by approximately 0.2137% or 24.1 m across all six methods. This is a finite-pressure truncation response within the input tables, not the physical missing-crust correction or P=0 radius.

### Stop rule and next physical gate

Every selected cutoff has P_cut>0. Even source row 1 is an internal isobar, not the vacuum surface. The R14 scan therefore does not admit a complete-star surface, radial \\(\\omega_0^2\\), or full-star Cowling spectrum. R13 synthetic tails remain reproducible historical diagnostics only. Continue to physical radial stability only after a sourced low-density crust EOS and a justified surface condition are available; alternatively, explicitly scope a finite-domain eigenproblem with a positive-pressure boundary and do not mislabel its roots as full-star f/p modes.


## 20. R15 full-cohort finite-pressure low-density cutoff diagnostic (2026-11-10)

R15 was executed as a separate EOS-preserving diagnostic over the frozen projected spectrum cohort. It measures how (M(r_{\rm cut})) and (r_{\rm cut}) vary as the integration endpoint is moved through exact existing native knots in the low-density prefix. It does not supply a physical outer-crust model or replace the mass-dependent f/p1 sequence.

### Frozen input and design

- 234 projected solver tables, 13 crust families, 3 admitted core families, 39 fixed crust-core pairs and 6 smooth matching prescriptions.
- Eight exact source-row cutoffs per table: `(1, 9, 17, 23, 52, 73, 116, 217)`.
- The operational segment is the contiguous prefix (n_b\le10^{-3}\,\mathrm{fm}^{-3}). All 234 projected tables have exactly 217 rows in this prefix and the same eight selected row indices.
- The historical `maxwell` token is still a smooth bridge; it is not a sharp density-jump Maxwell phase transition.

### Numerical completion and verification

Full-cohort native-density preflight Actions [42650145210](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42650145210) passed for 234/234 source tables and 1,872/1,872 selected knots. The full-cohort TOV Actions run [42650666976](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42650666976), exact source commit `16f44113464987c2043d1b6c70a028f936af16e8`, completed 1,872/1,872 integrations with all statuses `ok`.

Independent artifact revalidation confirmed:
- all 234 projected EOS hashes equal the frozen summary and are unchanged before/after;
- all 234 source-row-1 native-edge M/R anchors reproduce, with maximum relative errors (1.33\times10^{-15}) in mass and (1.22\times10^{-15}) in radius;
- all 312 pair×cutoff six-method summary cells are complete;
- all 1,872 profile SHA-256 and data-row counts match their manifest records;
- zero case failures and zero synthetic EOS tables.

The original 752-MB Actions artifact was repackaged without rerunning TOV in Actions run [42651779185](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42651779185). It yielded a core receipt/result artifact, four profile shards, and an index, each downloaded and independently verified.

### Descriptive outcome

At source row 217, relative to the same method/table source-row-1 finite-pressure isobar, the mean enclosed-mass change is (-0.00301608\%) and the mean isobar-radius change is (-3.934502\%). Core-family-conditioned mean radius responses are BBB_BHF2 (-3.475964\%), GM_GM1 (-4.622658\%), and RG_SLY4 (-3.704885\%). The maximum six-method spread at row 217 across complete pair cells is (0.00147311) percentage points in radius and (2.46122\times10^{-6}) percentage points in mass. These are descriptive contrasts in the fixed 1.3300106263-(M_\odot) reference cohort, not a universal causal response or a mass-sequence result.

All cutoffs have (P_{\rm cut}>0). The radius changes are internal-isobar coordinate changes, not complete-star radius or crust-thickness corrections; the enclosed-mass changes are not crust-mass estimates.

### Scientific next step and blocking gates

The 217→218 native density interval is sparse and crust-family-dependent. Across 13 representative crust families, (n_b(218)/n_b(217)=7.434)–(11.152). The representative two-point secant (Delta\\ln P/Delta\\ln n_b) is not a local sound speed, adiabatic index or causality test.

The full-cohort solver-native Hermite audit is now complete: Actions [42654562160](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42654562160) evaluated all 234 row-217→218 cells using the exact `cubic_Hermite_P_of_mu_raw_n_v1` construction. Every cell passes the Fortran-equivalent whole-cell (P''>0) and (\mu P''-P'\\ge0) checks; analytic `chi(t)` stationary-point enumeration finds zero interior roots, and `chi` increases monotonically over every cell. The global (c_s^2/c^2) range is (1.830688\\times10^{-4})–(2.950720\\times10^{-3}), with maximum endpoint First-Law closure residual (2.392\\times10^{-16}). Source hashes are unchanged; no EOS rows or TOV solves were added in this audit. Full evidence: `docs/R15_NATIVE_HERMITE_ROW217_218_AUDIT_20261110.md`.

R16 has now completed the full 234-table normalized-to-projected comparison; see `docs/R16_NORMALIZED_TO_PROJECTED_AUDIT_20261110.md`. The source q=n_B/mu_B sequence is already monotone in all 234 tables, so PAV does not pool knots. Pressure is instead reintegrated from dP/dmu_B=n_B. This exactly reproduces every frozen projected table, but the source pressure knots have 447 downward adjacent-secant transitions across all 234 tables. Thus a convex interpolant cannot preserve every source pressure knot exactly in the affected data sequences. Raw normalized Hermite cells reject 24,505/109,776 for stability/causal-margin conditions, while projected cells reject 0/109,776. This is a representation-level result, not a claim that the upstream microphysical EOS is physically unstable.


R15 does not lift MD-31 or MD-32:
- MD-31 still requires a sourced physical low-density crust EOS down to (P=0) and a validated vacuum surface.
- MD-32 remains a separately defined finite-positive-pressure eigenproblem only if its truncated-domain boundary condition is explicitly derived and justified.
- Direct radial (\omega_0^2) and complete-star mass-dependent Cowling f/p1 remain blocked. Barotropic-equilibrium (Gamma_1=\Gamma_{eq}) is distinct from frozen-composition (Gamma_1); the latter requires composition derivatives.

The authoritative R15 result record is `docs/R15_FULL_COHORT_TOV_RESULTS_20261110.md`; input preflight is `docs/R15_FULL_COHORT_NATIVE_DENSITY_PREFLIGHT_20261110.md`.

## 21. R17 constrained thermodynamic reconstruction comparison

R17 should treat the projection as a model choice with an explicit data-fidelity contract, not as a generic cleanup step.

### Baseline A: current q-preserving exact-integral construction

Preserve the normalized source (\mu_B) and (n_B) knots, use the nondecreasing (q=n_B/\mu_B), represent (q) linearly between knots, integrate (dP/d\mu_B=n_B) exactly, and reconstruct (\epsilon=\mu_B n_B-P). This is the current projected solver input and is the baseline to reproduce unchanged.

### Alternative reconstruction families to define

1. **Pressure-relaxed convex fit:** allow pressure knots to move while keeping source (\mu_B,n_B) knots fixed; enforce monotone/convex (P(\mu_B)), positive energy density and causal (c_s^2/c^2\le1). This family directly addresses the 447 source secant drops, but its endpoint pressure residual must be reported rather than hidden.
2. **Joint (P,n) reconstruction:** allow both pressure and density knots to move under (dP/d\mu_B=n_B), (P''>0), and causal-margin constraints. Explicitly report changes in (n_B), (P), (\epsilon), and (\mu_B); do not claim all source variables are preserved.
3. **Piecewise source segmentation:** where a phase/composition transition is supported by source metadata, test whether the offending knot sequence should be segmented at a justified interface rather than smoothed across it. A split is not licensed by a numerical residual alone.

### Objective and uncertainty rules

The current repository does not provide a complete statistical uncertainty model for all normalized EOS knots. Therefore:
- do not invent per-knot standard deviations or claim statistical optimality;
- expose objective weights and normalization scales as explicit sensitivity parameters;
- report a Pareto/error budget for pressure, density, energy, chemical potential, convexity, causality, and differentiability;
- keep exact byte identity of parent normalized/projected inputs separate from semantic table identity;
- store each candidate reconstruction as a separately versioned derived artifact with source-summary SHA, per-table parent hashes, algorithm version and full receipt.

### Acceptance and stopping conditions

- Baseline A must reproduce the current 234 projected CSVs at the existing (5\times10^{-12}) relative comparison threshold.
- Each alternative must pass per-cell (P'>0), (P''>0), (\mu P''-P'\ge0), (\epsilon>0), and the exact First-Law reconstruction.
- Any source pressure knot that cannot be preserved under convexity must appear in the explicit residual ledger; no silent knot replacement.
- Compare only source-backed low-density model/phase alternatives. Keep raw CompOSE package, normalized matched input, projected solver input and new derived candidates as distinct provenance layers.
- Do not rerun the full 1,872 TOV cases until at least two clearly defined, physically interpretable candidate reconstructions exist. Then select a small representative structural sensitivity slice before scaling to the 234-table cohort.

R16 did not run TOV or modes. MD-31 ((P=0) physical low-density EOS/surface) and MD-32 (finite-domain perturbation boundary) remain separate. Radial (\omega_0^2) and complete-star Cowling modes remain blocked.
