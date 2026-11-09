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
2. Run PYTHONPATH=python_frontend python -m pytest -q python_frontend/tests/test_mass_sequence.py.
3. Inspect the pilot artifact and manifest. Confirm the frozen 234-row input CSV hash and every projected EOS SHA-256.
4. Inspect each model's baseline reproduction status, branch status, target mass coverage, shooting residuals, W/V node counts and mass residuals.
5. Push pilot group is BBB_BHF2 × GPPVA_FSU2, six methods and eight requested masses. Validator requires a complete six-method 1.33-solar-mass paired cell and TOV M/R baseline reproduction.
6. If pilot passes, manually dispatch workflow with cohort=full. Full run computes all 234 rows with 48 log-Pc TOV samples per EOS. It may run up to six hours; do not claim completion until run and artifact are inspected.

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
