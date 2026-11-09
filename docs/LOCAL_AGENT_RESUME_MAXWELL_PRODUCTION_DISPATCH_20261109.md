# Local coding-agent resume prompt: native Maxwell production dispatch

Date: 2026-11-09
Repository: `cosmosapjw-quantum/NS_Cowling`
Feature branch: `research/maxwell-native-fortran-v1`
Draft PR: https://github.com/cosmosapjw-quantum/NS_Cowling/pull/3

## Authority and workspace safety

The exact code+workflow head exercised by CI is `22b0faa3c69150f3ae4fffcd2a6cae5d8231bb6c`. The exact code-bearing implementation commit is `f8678cb13a18f63bad9f8db65c4cd3e66e491ad3`.

The exact latest-head push and PR workflows both passed:

- Push run: https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42585358415; artifact `13658722733`; SHA-256 `26313f4315af74718f5fa9d1687658ef603c405807b2f7d44d27057fa38242a1`.
- PR run: https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42585366571; artifact `13659056245`; SHA-256 `12afce4e910c53b9a0dda180826bea2a98a5986562df0866bb7f9a26da219412`.

Before changing any local files:
1. Run `git status --short --branch`, `git rev-parse HEAD`, and `git branch --show-current`.
2. Preserve any local modifications. Do not reset, clean, checkout over, stash, or force-update user work.
3. Fetch `origin/research/maxwell-native-fortran-v1`. If the local branch is behind and clean, fast-forward only. If it diverges or is dirty, stop and report the exact divergence; do not reconcile by destructive reset.
4. Keep PR #3 open as draft. Do not merge it.
5. A pre-dispatch backup snapshot exists at `backup/maxwell-native-fortran-pre-production-dispatch-20261109`, pointing to the former implementation head `8d0aeb83a889d396794194338544d594c183e8b4`.

The initial implementation commit `f0723ecdfd90fcc5eda3ae27d73859f2b75e06f3` failed compile because the production generic wrapper omitted the optional smooth `eos_events(:)` argument. The fix is recorded in `f8678cb...`. Preserve that CI failure in the history; do not relabel it as a passing run.

## Research objective

Remove the previous smooth-EOS-only production boundary without flattening or interpolating through a sharp Maxwell density jump. Maintain a single type-resolved public entry point for smooth and one-interface Maxwell backgrounds, while preserving the existing low-level kernels and backward compatibility.

## Implemented public API

Source: `oscillation_source/solutions/src/production/ns_production_solvers.f90`.

```fortran
call solve_tov_star(eos, pc, options, smooth_profile)
call solve_tov_star(eos, pc, options, smooth_profile, eos_events)
call solve_tov_star(eos_high, eos_low, pc, pt, options, maxwell_star)

call solve_cowling_spectrum(smooth_background, omega_grid, mode_options, smooth_spectrum)
call solve_cowling_spectrum(smooth_profile, omega_grid, mode_options, smooth_spectrum)
call solve_cowling_spectrum(maxwell_star, omega_grid, mode_options, 'slow', maxwell_spectrum)
call solve_cowling_spectrum(maxwell_star, omega_grid, mode_options, 'rapid', maxwell_spectrum)
```

- Smooth TOV wraps `integrate_tov` and forwards the optional smooth pressure-stop list.
- Maxwell TOV wraps `integrate_maxwell_tov`, with distinct high/low EOS callbacks and exact `P_t` segmentation.
- Smooth Cowling spectrum wraps `scan_modes`; a `star_profile_t` can be copied into the existing `background_table_t` representation.
- Maxwell Cowling spectrum is implemented in `oscillation_source/solutions/src/modes/ns_maxwell_mode_spectrum.f90`. It takes one fixed conversion prescription per scan, scans an increasing positive angular-frequency grid, preserves trial statuses, builds brackets, uses bisection for sign/near brackets, reruns the shot at the root, and checks the final residual.
- A valid Maxwell scan means the requested finite-grid trials and any detected bracket candidates completed. It is not a completeness claim for the full spectrum.
- Root labels remain `unresolved`; no `f`, `p1`, or interface-mode identification is performed.

Production routes now use the generic API:
- Smooth TOV case adapter: `src/structure/ns_tov_cases.f90`.
- Normalized-table TOV adapter: `src/structure/ns_table_star_cases.f90`.
- Smooth mode cases: `src/modes/ns_mode_cases.f90`.
- Maxwell table case adapter: `src/structure/ns_maxwell_cases.f90`.
- Added CLI case `maxwell_table_spectrum`.

## Verified at exact code head f8678cb

Toolchain: GNU Fortran 16.2.0, FPM 0.13.0, flags:
`-std=f2018 -O0 -g -fcheck=all -ffree-line-length-none -fimplicit-none`.

Results:
- `cell_checks`: 12 passed.
- `maxwell_checks`: 22 passed, including smooth/Maxwell TOV dispatch, smooth profile Cowling dispatch, and slow/rapid Maxwell scan no-root path.
- `maxwell_table_checks`: 74 passed.
- Production `case_driver` build passed.
- Existing `maxwell_table_cowling` synthetic two-table CLI smoke passed.
- The analytic fixture small-jump comparison remains `dM/M=9.1343e-8`, `dR/R=2.8981e-8`.

The updated workflow now runs `maxwell_table_spectrum` on the generated analytic two-table fixture. It scans 1799, 1800 and 1801 Hz; residuals are approximately 0.989787, 0.993524 and 0.996295; all three trials succeed, with zero brackets and zero roots. Thus the production CLI no-root scan path is CI-verified, but the positive-root/bisection branch is not.

## Latest exact-head verification

At branch/workflow head `22b0faa3c69150f3ae4fffcd2a6cae5d8231bb6c`, runtime-check CI reports 12 cell, 22 analytic Maxwell and 74 normalized-table checks passed; production build and both Maxwell CLI cases passed. The synthetic scan is no-root only. The malformed workflow attempt at `42b5936...` and the earlier optional-argument compile failure at `f0723ec...` remain in history, followed by corrected commits and passing runs.

## Next research steps, in order

### 1. Verify the positive-root branch of the Maxwell scanner
The current fixture window around 1800 Hz has no residual sign change. Therefore `bisect` and the final root shot are implemented but not exercised by the current Maxwell spectrum test.

- Explore a broad enough angular-frequency grid for the existing analytic Maxwell fixture.
- Record every sample frequency, shot status and residual; do not infer a root from one small residual.
- Find an isolating sign bracket for slow and rapid conversion independently.
- Add assertions that the candidate is bracketed, the bisection root is valid, a fresh root shot is valid, and `abs(residual) <= residual_tol`.
- Add grid-refinement and root-tolerance refinement checks.
- If this fixture has no robust root, introduce a narrowly scoped deterministic seam for residual-bracket logic and test it against the real scanner interface. Do not invent or hardcode an eigenfrequency.
- Keep family label `unresolved` until phase-aware identity is established.

### 2. Independent Python/SciPy TOV cross-validation
Use the separate phase-aware Python/SciPy reference implementation, not a reimplementation of the Fortran equations in the test. Match EOS knots, SI constants, pressure cutoff, central pressure, and interface pressure. Compare (r_t), (m_t), surface (M), surface (R), and one-sided metric/thermodynamic values. Record tolerances and provenance.

### 3. Numerical convergence envelope
Vary ODE tolerances and max step size independently for:
- core TOV segment,
- envelope TOV segment,
- core Cowling segment,
- envelope Cowling segment,
- interface-transfer residual,
- root bisection tolerance once a positive bracket fixture exists.

Do not use a single small residual as a convergence study.

### 4. Admission negative tests
Verify rejection and precise status for:
- (P_t) outside either phase-table common domain,
- mismatched (mu_B(P_t)),
- reversed/nonpositive (Delta n_B) or (Deltaepsilon),
- first-law jump-closure failure,
- failed phase callback during integration,
- incomplete frequency scan with failed trials,
- zero/duplicate/unsorted/nonfinite frequencies,
- unresolved near-zero sample without an isolating bracket.

### 5. Physical data gate
Only then attach a realistic source-backed CompOSE or other hybrid EOS pair. Preserve source file, release/version, units, checksums, transformations, phase labels, crust/core matching and the declared coexistence pressure. Do not call the synthetic tables realistic data.

### 6. Later branches
Central-pressure continuation and stable M-R branches, radial stability, twin-star analysis, phase-aware mode identity/continuation and observational inference remain downstream. The Seidov ratio is only a local onset diagnostic.

## Commands for local verification

From `oscillation_source/solutions`:

```bash
fpm test --compiler gfortran --flag '-std=f2018 -O0 -g -fcheck=all -ffree-line-length-none -fimplicit-none'
fpm build --compiler gfortran --flag '-std=f2018 -O0 -g -fcheck=all -ffree-line-length-none -fimplicit-none'
```

Then exercise both generated fixture cases:
- `case_driver maxwell_table_cowling ../../dist/maxwell_table_fixture.nml`
- `case_driver maxwell_table_spectrum ../../dist/maxwell_table_fixture.nml`

The fixture is generated by `maxwell_table_checks`. Assert `result_valid` and `cowling_shot_valid` for the first case; for the spectrum case assert `spectrum_scan_valid`, `spectrum_all_trials_succeeded`, and `spectrum_valid`. Also inspect the trial records and candidate/root fields. A true `spectrum_valid` without a root is a successful finite-grid no-root scan, not an eigenmode result.

## Non-claims to preserve

- No source-backed realistic CompOSE hybrid EOS has been executed in this implementation.
- The current synthetic table fixture is an integration regression fixture only.
- The Maxwell root-refinement branch needs a positive-bracket test.
- No root-family labels, f/p1 identity, interface-mode classification or spectrum completeness is established.
- No maximum mass, stable branch, twin-star, radial-stability sequence or astrophysical inference.
- No full-GR metric perturbations, rotation, magnetic field, elasticity, superfluidity, composition stratification or GW damping.
- The 234-row v5 smooth-matching results and old smooth `maxwell` bridge remain physically separate from sharp density-jump results.
- The separate 33-path v5 overlay remains unapplied unless independently proven otherwise on the exact repository tree.

## Completion protocol

After each independently passing unit:
1. record the exact Git commit and changed file list;
2. run the focused test and full relevant FPM/CLI verification;
3. push a fast-forward commit to `research/maxwell-native-fortran-v1`;
4. retain exact-head Actions logs and artifact digest;
5. move/advance a named backup branch only after recording its target SHA;
6. update this resume prompt, the DAG, the implementation guide and adversarial audit;
7. update PR #3 metadata without merging or marking it ready for review.


## Superseding authority and positive-root regression (2026-11-09)

This section supersedes earlier authority/result statements above.

- Exact code/test commit: `9e92655a8277922f312f4c31df8e331e8c12df1b`.
- Exact strict workflow head tested by both push and PR: `ec4527f0d518f6ba19f57eb42d8a8eb00b98dfca`.
- Push [42586438396](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42586438396), artifact 13658709148, SHA-256 `02d1673cd82045aaacbae8125c74c7ce7b8896b7abe6def70e81df07331c3d04`.
- PR [42586446195](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42586446195), artifact 13659137810, SHA-256 `478fb960cf8a993bb92bf2bbc544f17160b854e0f0f5ac5a28cc0b14a9ccaeaf`.

The latest CI requires exactly 2 valid roots for each conversion prescription on a 41-point 100-20,000 Hz grid. Slow roots: 1764.8989380477 Hz (residual -6.9575e-7), 4937.6449322607 Hz (9.3334e-8). Rapid roots: 4937.0162338857 Hz (8.5331e-9), 19493.6326018255 Hz (3.8951e-7). All mode-family identities remain unresolved; these are synthetic analytic-table results, not realistic star predictions.

Next work: perform independent root-frequency convergence by varying ODE tolerance, maximum step, bisection tolerance and grid density separately, then compare the sharp-interface Python/SciPy reference with identical table EOS and background. Do not assign mode identity by root ordering.


## Superseding verified state: root convergence and Python/SciPy reference (2026-11-09)

This section supersedes earlier status/next-step statements above.

- Exact Fortran/table fixture commit: `9d3ef4bfe8b018b17887504d60e25d1c94971b6e`.
- Exact validator head: `afcdf948f0ca6844a6c3402f66e184a46aa78978`.
- Push [42588162457](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42588162457): artifact 13659334940, SHA-256 `2f394e17b8f2afa2618322c3c634f63fb6faff81090fe14c700d66f374f165d5`.
- PR [42588170532](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42588170532): artifact 13660500740, SHA-256 `05d6786c5a5e82e177646a7db72efc50d485d54d0fac45ae1cec8a499666fd8e`.

Strict CI results: 12 cell, 22 analytic Maxwell, 74 normalized-table assertions; production build; baseline slow/rapid positive-root spectra; ODE-refined slow/rapid spectra; denser-grid/tighter-bisection slow/rapid spectra; independent Python/SciPy reference. The strict validator gates are:
- ODE refinement root shift <= 1e-4; measured maximum 3.2810e-5.
- Frequency-grid/root-tolerance refinement shift <= 1e-7; measured maximum 2.6561e-10.
- Python TOV relative differences <= 1e-6; measured M 5.4709e-10, R 2.1517e-10, r_t 2.6674e-10.
- Python Cowling root-set difference <= 1e-4; measured maximum 5.2046e-5.

Reference frequency set, sorted only for numerical comparison:
- Rapid native 4937.0162338857 / Python 4937.0455932778 Hz; native 19493.6326018255 / Python 19494.6472160782 Hz.
- Slow native 1764.8989380477 / Python 1764.9225737103 Hz; native 4937.6449322607 / Python 4937.6746742499 Hz.

No mode identities are assigned. The root list is not a complete-spectrum claim and the fixture is synthetic. A source-backed realistic CompOSE hybrid EOS remains absent.

## Next local research tasks

1. Implement phase-aware mode identity using segment-aware node counting, a consistent normalization/overlap, and continuation in frequency. Treat the slow/rapid spectra as separate boundary-value problems; do not pair roots by order.
2. Add negative admission tests for out-of-range transition pressure, mismatched chemical potential, reversed jumps, first-law closure failure, callback failure, failed scan trial and non-bracketed near-zero sample.
3. Compare one-sided Cowling profiles and interface-transfer diagnostics pointwise against Python/SciPy at matched frequencies and tolerances.
4. Construct a source-backed CompOSE hybrid EOS pair with provenance, units, phase labels, crust/core matching, coexistence pressure and checksum; only then run the production solvers.
5. Build central-pressure sequences, M-R branches, turning-point/radial stability and twin-star diagnostics as separate gates.

Do not merge PR #3 or remove the draft label. Preserve old failed Actions attempts and the pre-dispatch backup branches as provenance.
