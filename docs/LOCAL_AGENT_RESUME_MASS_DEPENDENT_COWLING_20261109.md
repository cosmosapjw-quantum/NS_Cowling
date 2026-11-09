# Local coding-agent resume: mass-dependent smooth-matching Cowling spectra

## Authority and safety

- Repository: cosmosapjw-quantum/NS_Cowling
- Branch: research/mass-dependent-cowling-v1
- Branch parent at creation: a09ec30d74b363652d2d1aefe990557c468679dd
- Fetch and report the live branch HEAD before every mutation; never assume this prompt remains current.
- PR #3 is open/draft/unmerged on research/maxwell-native-fortran-v1. Do not merge, close, retarget or mark it ready.
- PR #2 is also open/draft/unmerged. Do not modify it.
- No user workstation checkout was edited. Remote Desktop Commander reported the desktop device offline.
- Preserve all historical Maxwell backup branches and failed Actions attempts.

## Scientific design

This is the historical smooth crust-core matching cohort, not the sharp first-order Maxwell density-jump track. It contains 13 crusts × 3 admitted cores × 6 smooth methods = 234 method rows and 39 fixed core/crust pairs. The old spectrum file is only at 1.3300106 solar masses; the reuse audit reports recomputed_modes=0.

Mass targets: 0.8, 1.0, 1.2, 1.33, 1.4, 1.6, 1.8 and 2.0 solar masses.

The method token maxwell means a smooth near-pressure-plateau bridge with no latent heat or density jump.

## New implementation files

- python_frontend/compose_cowling/mass_sequence.py: stable-branch extraction, hard barriers at failed TOV samples, adjacent-Pc mass brackets, complete six-method paired summaries.
- python_frontend/tests/test_mass_sequence.py: branch maximum, gap barrier, log interpolation and no-imputation tests.
- python_frontend/compose_cowling/solver.py: optional explicit pc_bracket_pa. Old default is retained.
- scripts/run_mass_dependent_spectrum.py: TOV scan, stable branch, target-mass solve, Cowling spectra, node checks and data/provenance manifests.
- scripts/validate_mass_dependent_spectrum.py: pilot/full gates.
- .github/workflows/mass-dependent-cowling.yml: push-triggered six-model pilot and manual pilot/full dispatch.
- docs/MASS_DEPENDENT_COWLING_RESEARCH_PLAN_20261109.md
- docs/MASS_DEPENDENT_COWLING_DAG_20261109.md

## First required checks

1. Fetch current branch HEAD and the workflow run tied to that SHA.
2. Run PYTHONPATH=python_frontend python -m pytest -q python_frontend/tests/test_mass_sequence.py and compile the native `cell_solve` executable.
3. Inspect the pilot artifact and manifest. Confirm the frozen 234-row input CSV hash and every projected EOS SHA-256. The failed first attempt is Actions run 42598463512 / artifact 13664393404; it produced zero valid spectra and is a failed integration test, not scientific evidence.
4. Inspect each model's baseline reproduction status, branch status, target mass coverage, shooting residuals, W/V node counts and mass residuals. The first run's M/R mismatch was traced to using `case_driver real_table_star` (r0=1 m, rtol=1e-8) against v5 `cell_solve` metadata (r0=10 m, rtol=1e-9, step_cap=0.05). The new runner must use `solve_cowling_tov_star` / `cell_solve task=star` and identical v5 settings for N1.
5. Push pilot group is BBB_BHF2 × GPPVA_FSU2, six methods and eight requested masses. `cell_solve` has a configurable `r_limit_m`; its old 1e6-m default is preserved, while the exploratory scan uses 1e8 m. Confirm the radius-limit setting is carried through mass-root and mode matching calls. Validator requires a complete six-method 1.33-solar-mass paired cell and TOV M/R baseline reproduction.
6. Do not dispatch full cohort until a new pilot passes. If pilot passes, manually dispatch workflow with cohort=full. Full run computes all 234 rows with 48 log-Pc TOV samples per EOS. It may run up to six hours; do not claim completion until run and artifact are inspected.

## Required analysis after full execution

- Build a 39-pair × mass coverage table. Distinguish first_turning_point_resolved, rising_branch_truncated and sequence_gap_truncated.
- For y in {f,p1}, calculate
  \[
  \Delta_m^{(y)}(M;c,s)=100[y_m(M;c,s)/\bar y_6(M;c,s)-1].
  \]
- At every mass report complete-pair denominator, median, 90th percentile and maximum six-method range.
- Provide core/crust-resolved maps, method sign reversals, and mass of peak p1 sensitivity.
- Test, do not assume, whether core family shifts the peak sensitivity mass or method-response shape.
- Compare low-mass trend with Canullán-Pascual et al. (2025); merely reproducing low-mass p1 amplification is not the novelty.
- Compare representative Cowling and full-GR literature results. Keep Cowling systematic offsets separate from paired matching shifts.
- Add direct published matching baselines: simple concatenation only as a deliberately non-thermodynamic diagnostic control; linear P(epsilon); Canullán-Pascual et al. (2025) thermodynamically consistent P(epsilon) interpolation; and CUTER-like or otherwise documented thermodynamic crust reconstruction. Preserve provenance and run first-law, chemical-potential, causality and endpoint gates. Do not mislabel an approximate reconstruction as CUTER.
- Refine low-mass tail, reference mass and high-mass cases in TOV/Cowling tolerances, maximum step and root brackets. Require frequency movement <=0.1 of the smallest reported method effect.
- Only after these gates add interface radius/enclosed mass, crust thickness/mass, tidal Love number/Lambda, source-backed CompOSE EOS and observational inference.

## Claim boundaries

- The 234 rows form a finite crossed design, not 234 independent EOS samples or a posterior.
- maxwell is a smooth bridge, not a sharp density jump.
- The current calculation is ell=2 relativistic Cowling with a finite positive pressure cutoff, not full-GR complex QNMs or damping times.
- The first M(Pc) turning point is a barotropic branch/stability proxy; a radial pulsation calculation is required before an unqualified radial-stability claim.
- Current output includes M,R and f,p1. It does not yet include crust mass/thickness, tidal Lambda or a sensitivity kernel.
- Failed/out-of-branch targets remain unavailable; no interpolation or imputation.

## Safe stopping rule

If pilot CI fails, identify the exact failure class: EOS pressure domain, TOV integration, sequence gap, branch resolution, target-mass bracket, Cowling frequency bracket, node identity, residual or baseline mismatch. Do not widen frequency brackets silently or loosen acceptance tolerances to force green. Preserve logs/artifacts and repair the smallest upstream cause.

## Publication gate

JKPS/EPJA is a plausible first target only if full mass-dependent paired response, explicit coverage, numerical refinement, and independent benchmarks reveal a robust core/crust-conditioned result beyond already published low-mass p1 sensitivity. A&A/MNRAS additionally requires realistic EOS provenance and observational/tidal inference. PRD requires a new relativistic perturbation/kernel/full-GR contribution.


## Source identity and exact next gate

The v5 source archive checksum recorded in release_summary.json is b2c0dcaf6a5b07956cf8d33c39651d4ade32f44b4c49b33aa08fd84fb458185b. It is a reference receipt, not an independently verified archive checksum in this session. The runtime v5 backend hashes the current 35 source files from the branch, their relative paths, and compile_manifest.json. It records a backend build key and manifest SHA in each model record. Byte identity against the archived ZIP/TAR remains UNVERIFIED.

After push:
1. inspect the exact-head Actions run and confirm it compiles calculation/methods/support, not solutions;
2. check the runtime v5 backend build key and manifest SHA in the artifact;
3. compare frozen M/R at each stored Pc before considering mode results;
4. only then inspect the first fold candidate and complete pair coverage.


## Current continuation state: branch topology is the active blocker (2026-11-09)

Latest fully tested code/theory head before the candidate-segment increment: d1d01e6b114d97cd496640037aafaf7ff7be126e. Backup refs:
- backup/mass-dependent-pre-cowling-variation-20261109 -> ef738a620c440fea7e91d4c10da13d2f2e30c93a
- backup/mass-dependent-pre-segment-audit-20261109 -> d1d01e6b114d97cd496640037aafaf7ff7be126e

The exact-head Actions run 42610037961 (job 129135805050; artifact 13668724637; zip SHA-256 fb62570658e44f9847d4a08f73174e36f14c091b6960b3b847206232a73b1534) completed with failure at the physical pilot validator. The two helper/mock suites passed (6 tests each). The validator reports: required mass 1.0 lacks a complete six-method pair. All six models reproduce the frozen v5 TOV M/R at the reference Pc to near machine precision, but the 24-point Pc scan's first rising-proxy maximum is only 0.7432-0.7434 Msun, so every requested mass >=0.8 Msun is outside the initial rising proxy. Valid f/p1 roots = 0. No scientific mass-dependent spectrum is admitted.

The 24-point scan has a second, higher-density rising interval ending at about 1.915 Msun, after an extreme excursion in M and R. This interval is not labeled stable or used for mode roots. The new implementation adds extract_candidate_mass_segments, candidate_mass_segments.csv, explicit radial_stability_status=unresolved_requires_radial_mode_spectrum, and a 96-point six-model pilot. It also renames branch output/diagnostics to avoid calling the first local peak Mmax. Inspect the exact new HEAD and its Actions before trusting the implementation.

## Immediate agent tasks

1. Fetch the live branch ref and exact-head Actions. Inspect code, branch-proxy and test changes; do not assume the current handoff's SHA is still head.
2. Confirm the four candidate-segment unit tests: rising/falling/rising split; shared maximum plateau; hard barriers at failed or omitted indices; unresolved radial stability labels.
3. Inspect candidate_mass_segments.csv and tov_sequence.csv from the 96-point artifact for all six methods. Report local maxima/minima, M/R/Pc, failed integration points, and grid coverage. Do not label any candidate segment radially stable.
4. Require 96-to-192 resolution convergence for the first peak and rapid M/R excursion. Do not infer physical EOS discontinuity from a coarse grid.
5. Reproduce selected points independently with exactly the native raw-knot P(mu) cubic-Hermite EOS interpolation. A generic P(epsilon) spline is not an independent comparator for this code path.
6. Implement or integrate a relativistic radial pulsation Sturm-Liouville solver. In geometric units the standard two-variable equations for xi and DeltaP, centre BC, and zero-pressure surface BC are given in docs/MASS_SEQUENCE_BRANCH_STABILITY_AUDIT_20261109.md. The v5 finite-positive-pressure endpoint needs either a low-density EOS continuation to P=0 or a separately derived finite-cutoff boundary condition.
7. Declare the perturbation closure. Gamma_eq follows the equilibrium barotrope; Gamma_1=(epsilon+P)/P*(partial P/partial epsilon)_(s,Y_i) is the frozen-composition adiabatic response. If composition derivatives are unavailable, only label results as barotropic-equilibrium radial stability, not general physical stability.
8. Track the nodeless omega_0^2 and its zero crossings along Pc. Only after the radial eigenvalue/node gate passes should mass-target f/p1 spectra be attempted on each radially stable interval.
9. For a true sharp Maxwell interface, implement slow/rapid/reactive radial junction conditions separately from the historical smooth six-method cohort.
10. Do not dispatch the full 234-EOS cohort until the six-model high-resolution pilot passes its branch, radial-stability, M/R, mass-residual, shooting and node gates.

## Scientific reference boundary

- Alford, Harris & Sachdeva 2017, arXiv:1705.09880: direct radial Sturm-Liouville eigenvalues are the fundamental stability test; BTM turning-point heuristic can differ for exotic low-compactness branches.
- Essick 2024, arXiv:2405.05395: EOS-dependent stability recovery/secondary branches can occur in selected low-pressure, low-compactness sequences.
- Rau & Salaben, arXiv:2309.08540v2: distinguishes equilibrium Gamma from frozen-composition Gamma_1 and derives interface junction conditions for genuine first-order transitions.
- Kokkotas & Ruoff 2001, DOI 10.1051/0004-6361:20000216: published radial-mode framework and reference frequencies.

The theory/mock commit d1d01e6 and this candidate-segment increment are separate from PR #3. PR #3 must remain open/draft/unmerged. No local workstation checkout was modified. Dropbox artifacts were inspected read-only; no Dropbox mutation was performed.


## Exact 96-point artifact readback and immediate 192-point gate (2026-11-09)

Current tested candidate-segment implementation HEAD: dd5a33e521ecde554945c82b79cfd534ed91b306. Backup before 192-point refinement: backup/mass-dependent-96point-pilot-20261109 -> dd5a33e521ecde554945c82b79cfd534ed91b306.

Exact Actions run 42611801783 failed only at the physical validator after all six TOV rows and both test suites completed. Artifact 13670407745, ZIP SHA-256 ecafc110e8b9d0196dfb4bb341da2ee54273f5d61767cca7c36e502a40737b7d. All six EOS/method rows completed 96/96 TOV points with no failed samples, and the frozen baseline M/R anchor reproduces to near machine precision. Unit tests: 10 mass-sequence tests and 6 fixed-mass/Cowling mock tests passed. The 48 requested spectra are unavailable; complete paired cells=0; required 1.0 Msun validation fails correctly.

The 24-point first maximum was near 0.7433 Msun / 0.172 MeV fm^-3 / 1233 km; the 96-point first maximum is 0.76502-0.76511 Msun / 0.22676-0.22733 MeV fm^-3 / 1879.8-1889 km. Thus 24-to-96 is not converged. Five methods have four monotonic pieces; the smooth maxwell bridge has six, including an extra max/min at Pc about 1.116/1.898 MeV fm^-3, M about 0.21684/0.20288 Msun and R about 17.86/17.31 km. The bridge has no density jump; its metadata gives cs2_min about 9.07e-4. No piece is called stable.

The workflow is changed to 192 Pc samples on the six-model pilot. After push, inspect the exact-head run/artifact and compare first turning-point M/R/Pc plus every candidate fold. Do not use the secondary/higher-density intervals for Cowling roots until direct radial fundamental-mode omega_0^2, an explicit Gamma_1 closure and a valid zero-pressure/finite-cutoff surface condition have been validated.


## Exact 192-point artifact readback and next action (2026-11-09)

Current 192-point pilot commit: 40177800d4484a19149485581e0f2d5bd69594a9. Backup: backup/mass-dependent-192point-topology-20261109 -> 40177800d4484a19149485581e0f2d5bd69594a9.

Actions 42615970377 completed the six-model 192-point run. Artifact 13673520039, ZIP SHA-256 e460c295af38a70444e6139fb82aa037b3be4718e29d172db7be31d358abbff5. All 1152 TOV integrations were valid; the six reference M/R anchors reproduced to machine precision. The 10 branch helper tests and 6 fixed-mass/Cowling mock tests passed. The physical validator failed correctly because all 0/48 spectra are available and the 1.0 Msun six-method cell is absent.

The first sampled max is M=0.7652434-0.7652747 Msun, Pc=0.2248758-0.2254394 MeV fm^-3, R=1849.995-1858.853 km. The 96-to-192 mass shift is small, but R shifts about 1.6% and the first subsequent minimum R in css shifts about 418 to 195 km. The maxwell-labelled smooth bridge retains an extra maximum/minimum pair (M about 0.21795/0.20282 Msun at Pc about 1.0415/1.9619 MeV fm^-3); it is not an exact density-jump transition and is not a stability label.

The next workflow is configured for 384 central-pressure points and adds a SciPy DOP853 TOV reference using the exact same raw-knot cubic-Hermite P(mu) EOS. Five selected Pc values compare M/R against native v5 cell_solve with a predeclared 1e-6 relative criterion. Inspect both the 384-point topology and SciPy reference artifact before promoting any structural claim. Radial omega_0^2, closure Gamma_1, and surface condition remain unimplemented; no mass-dependent Cowling mode result is admitted.


## Exact 384-point result and current SciPy repair (2026-11-09)

Current branch before the repair commit: 6fb68cf1204b9656aa526f437a380e3c2e340a58. Backup created: backup/mass-dependent-384point-scipy-diagnostic-20261109 -> 6fb68cf1204b9656aa526f437a380e3c2e340a58.

Actions 42616463263 completed:
- 384/384 TOV points for each of six EOS/method rows (2304 total), no failed integrations;
- all frozen v5 baseline M/R reproductions pass;
- 10 mass-sequence tests, 6 fixed-mass/Cowling mock tests and 4 Hermite-EOS tests pass;
- no admitted f/p1 spectra: 0/48, no complete six-method target-mass cell.
Artifact 13671829050, SHA-256 a8bd61678ec74e1dea5451718f11386b44bb40162cad6ab29e10f1fe6268befc.

First sampled M(Pc) maximum:
- M = 0.7652866–0.7652926 Msun;
- Pc = 0.2239479–0.2245094 MeV fm^-3;
- R = 1835.56–1844.27 km.
First maximum mass differs by only about 2.6e-5 relative from the 192-point estimate, but R moves about 0.75%. The first subsequent minimum remains radius-unconverged: 384-point R spans 258–364 km across methods and differs non-monotonically from 192-point results.

The maxwell compatibility smooth bridge has an extra M(Pc) max/min pair:
- max M=0.2180954 Msun, Pc=1.0606304 MeV fm^-3, R=17.8517 km;
- min M=0.2028434 Msun, Pc=1.994654 MeV fm^-3, R=17.1793 km.
It has zero density jump and is not an exact Maxwell construction.

Initial independent SciPy TOV comparison run:
- 4 EOS-Hermite tests passed;
- five selected direct-r DOP853 vs native Fortran cases;
- maximum relative mass difference 3.60164e-7, radius difference 1.46572e-6;
- radius criterion 1e-6 failed narrowly. This is a numerical comparison concern, not a physical refutation; no tolerance was widened.

The repair now underway writes a diagnostic receipt before acceptance and records all per-case differences for:
1. coarse SciPy vs default native;
2. tight SciPy vs default native;
3. tight SciPy vs tight native;
4. native baseline vs native refinement.
SciPy rtol/max step and native q-step are refined. The 1e-6 acceptance remains unchanged. Inspect the exact-head CI artifact; do not call independent TOV validation PASS unless the tight comparison passes.

If the refined SciPy/native comparison passes, close the TOV structural cross-check only. If it still fails, preserve per-case receipts and classify the mismatch (surface event, EOS state, centre expansion, ODE integration) before any further repair. Do not widen the criterion. The radial eigenvalue solver remains unimplemented, so no branch is radially stable and no new mass-dependent f/p1 prediction is admitted.


## Continuation R10: exact-head receipt readback and CLI schema defect (2026-11-09)

Live branch before this repair: `57511496a57e1a1cee6bafc5254de63f4f5ab514`. Exact-head Actions run: `42617251360`; artifact: `13673098470`; ZIP SHA-256: `0d376752ba610819d3ca59fceedbf50ea5eebc5d03b3f440e3b87135546b3964`.

Readback:
- The 384-point six-model run completed 2304/2304 native TOV integrations with no failed samples.
- All six frozen-v5 reference M/R anchors reproduced.
- 10 mass-sequence tests, 6 fixed-mass/Cowling mock tests and 4 Hermite-EOS tests passed.
- Independent SciPy receipt status: `PASS_REFINED_SAME_EOS_INDEPENDENT_TOV`.
- Maximum tight SciPy/native relative differences: mass `1.6232670240598032e-8`, radius `1.1432510017123261e-7`, both below the unchanged `1e-6` threshold.
- Maximum native refinement differences: mass `1.299910193106295e-10`, radius `7.057766721929021e-11`.
- All this is a structural TOV comparison only. No radial stability or Cowling mode claim follows.

The workflow status was failure because the CLI report raised `KeyError: 'max_relative_mass_difference'` after writing the passing v2 receipt. The stored schema uses nested keys `max_relative_difference_scipy_tight_vs_native_tight`, `max_relative_difference_coarse_vs_native_reference`, and `max_relative_native_refinement`. The repair adds a `cli_summary` adapter and a regression test for `main()`; no acceptance tolerance is changed.

### Next exact-head checks

1. Inspect the live branch HEAD and PR #3 status before any mutation.
2. Confirm the CLI regression test passes and the summary prints all three receipt groups.
3. Inspect the new exact-head Actions result. The independent TOV step should complete; the subsequent physical mass-pair validator may still hold/fail correctly because 0/48 f/p1 roots and no complete 1.0 Msun six-method cell exist. Do not call the entire pipeline green unless the run says success.
4. Preserve the first rising-proxy endpoint near 0.7653 Msun as the target-coverage limit for that candidate segment. Do not extrapolate to 1.0 Msun and do not infer radial stability for the later high-density rising segment.
5. The next physics implementation is the direct radial pulsation eigenproblem for the nodeless fundamental omega_0^2. Explicitly choose the equilibrium Gamma versus frozen-composition Gamma_1 closure, and resolve the finite-positive-pressure surface condition. Without composition derivatives, make only the barotropic-equilibrium stability claim.
6. Only after the radial eigenvalue/node and boundary-condition gates pass should Cowling f/p1 roots be attempted on the newly classified stable intervals. Keep the true sharp-Maxwell slow/rapid junction laws separate from the six historical smooth bridges.

Dropbox backup target is the existing `/TOV_Cowling_Solver_Backups` folder. A dated checkpoint document should record the commit, Actions run, artifact digest, numerical receipt, and remaining science hold. PR #3 remains open/draft/unmerged and must not be altered.
