# C20 — native TOV PCHIP and batched Cowling geometry

Date: 2026-11-23 (Asia/Seoul).
Baseline: C19 commit 122c502d1ef5e55fd890f47d86444df725987032.
Branch: research/cowling-c20-native-tov-pchip-20261123.

## Decision and verified result

This scoped migration is implementation-verified on the original M2/M3 campaign.
TOV interpolation coefficients, dense radius/mass queries, and construction of
the five Cowling background fields now run in Fortran on the native source path.
Python packs source data, selects physical phases, invokes native functions, and
serializes evidence. Existing source-contract EOS state evaluation remains in
Python and is a separate remaining migration.

The complete scientific output for C19 serial, C20 serial and C20 two-process
runs is byte-identical: M2.json, M3.json, and both compressed raw scan payloads.
This includes40 stars,4,440 frequency shots,287 root candidates and80 f/p1 decisions.
No rejection was silently promoted or dropped.

| Configuration | Wall time(s) | Speed ratio |
|---|---:|---:|
| C19,1 process, fresh baseline |65.015313|1.000|
| C20,1 process |48.677465|1.3356x vs C19 serial|
| C20,2 isolated processes |24.820273|1.9612x vs C20 serial;2.6194x vs C19 serial|

Serial wall time fell 25.13%. Each number is one full measured campaign on the
same host, with the shared library built before timing; no statistical timing
interval or portable speed guarantee is claimed. The slower second model is not
parallelized by DVODE threads: each model's 20 densities run in its own OS process.

The separately timed Cowling coefficient construction fell from
29.177644303 s to 13.026055058 s (M2+M3 serial sum).
TOV integration and original source EOS sampling remain substantial work.
See data/cowling_c20_native_tov_pchip_20261123/PERFORMANCE_FINAL.json.

## Scientific invariants and identity

Original immutable archive:
data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip

SHA256:
299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04

Models: SPG_M2 Maxwell/rapid and SPG_M3 CSS/slow.
Closure remains cold equilibrium Gamma1=Gamma_eq. Signature(-,+,+,+).
Geometrized radii/masses are km; P,epsilon supplied to TOV are km^-2.
Source EOS data retain MeV/fm^3; no sound-speed derivative is interpolated
across a physical Maxwell jump.

C17/C19 controls remain:
TOV rtol3e-10,atol3e-11,max_step_h.004,native_samples2049,surface_h3e-8;
Cowling rtol1e-9,atol1e-12,max_step_h.002,root_xtol_hz1e-5;
4097geometry samples,360profile samples; source restart jump.01.
Normalized determinant and independent full-state join gate remain1e-5.
W/V node identity,12% frequency step, high/low phase overlap.95 remain unchanged.

| Model | Admitted f branches | Admitted p1 branches |
|---|---:|---:|
| M2 |20/20|16/20|
| M3 |16/20|7/20|

M3 factor 1.05 p1 high-phase overlap rejection remains intact.
Predeclared migration comparators were abs frequency<=.005Hz (C13) and
relative background<=2e-7 (C15). The actual final differences are0.0 for
all compared frequencies, backgrounds, overlaps and scan determinants.
The stronger complete-file comparison also passed.

The scientific manifest itself has SHA256:
dd033582b765c31828165d5dd9c83670de256f1610d539f765e87a1409e90254

Payload SHA256:
- M2.json: a3be0502ce6bdcb496e986882946254837b0faf01d18c650051f0b099927bef4
- M3.json: 21776b03975ef07f13d9423017ddd56eca0db7485e53f36b95a3e863f3adb3bd
- M2_all_scan_points.json.gz: 265e7359e6a09944f24a0c7a6e66f85d94e92d5303097bc61052fbbd6aafd323
- M3_all_scan_points.json.gz: 8ce744dc3f1afdace1e961f4e590246d129c6369abc22a8d670e97f7d3bc833d

## Formula-to-code map

| Mathematical operation | Fortran implementation |
|---|---|
| Shape-preserving slopes with weighted harmonic mean and one-sided limiter | compose_eos_cowling_tov_geometry:tov_native_pchip_slopes |
| Actual derivative at final knot, matching original SciPy query semantics | same routine, at_knots=1 |
| Radius/mass dense output at enthalpy queries | tov_native_phase_values |
| Phi=.5log(1-2M/R)-h+h_surface | tov_native_geometry |
| Lambda=-.5log(1-2m/r) | tov_native_geometry |
| Phi'=(m+4pi*r^3*P)/(r(r-2m)) | tov_native_geometry |
| Existing five-field geometry PCHIP | legacy API in compose_eos_cowling_pchip |

The input phase is selected before querying. Native batch construction retains
exact h endpoints and the existing high-side geometry choice at a shared
interface. The source EOS side remains selected independently and explicitly.

PCHIP follows the installed SciPy 1.17.0 definitions:
https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.interpolate.PchipInterpolator.html

The exact floating-point evaluation order was checked against installed
_cubic.py and tagged source:
https://github.com/scipy/scipy/blob/v1.17.0/scipy/interpolate/_ppoly.pyx
Git blob SHA5bcc2b911a0b11f5bb387e88473acdc5fb278e98.

The TOV path preserves reciprocal-of-normalized-harmonic-mean arithmetic,
CubicHermiteSpline coefficient construction, and ascending-power polynomial
summation. It does not substitute Horner evaluation for the historical path.
Native log1p retains weak-field accuracy. Native libm pow(r,3.0) preserves the
existing Python float cubing convention. These are Fortran calls into the
platform mathematics library; no Python interpolation callback is introduced.

## Why the first migration was rejected

The first mathematically equivalent Hermite implementation preserved all 80
selected f/p1 branch decisions but changed 44 of 287 raw candidate admissions.
Matched accepted frequencies differed by at most 0.000416572 Hz and backgrounds
by less than 4e-8 relative, yet this did not preserve the complete source record.
The first timing 45.805375 s is therefore not used as the final speed claim.

A bounded diagnostic isolated:
1. harmonic-mean rounding;
2. raw final slope versus evaluated final-knot derivative;
3. Hermite versus ascending-power dense evaluation;
4. integer-multiply versus Python-float radius cubing;
5. high-side versus low-side endpoint geometry within roundoff.

Three affected source cases(M2 factor .55, M3 factor .495, M3 factor .75) were replayed.
After repair, their complete TOV arrays, dense queries, geometry, geometry
slopes and candidate rows matched C19 exactly. The final full campaign then
matched all four science files exactly. No tolerance or root gate was relaxed.

The rejected runs, comparison failure and original shared binary are preserved
in C20_FIRST_REPLAY_REJECTED.tar.gz in the backup directory. Its retained binary
SHA256 is76dd38e4f969eaeb72962e3db0b5797d9159828176675a28facf6d3f7e93fe8c.
The final report does not claim that these preliminary results are certified.

## Memory and process behavior

Native phase objects store grid(n),values(2,n),slopes(2,n):5n doubles.
The old phase plus SciPy polynomial stores2n+8(n-1)+n=11n-8 doubles.
The new evaluator forms local coefficients with O(1) temporary workspace.

Actual source phase lengths gave:
- M2 factor .55:888760->404040bytes;
- M3 factor .495:359440->163440bytes;
- M3 factor .75:526288->239280bytes.

That is about 54.5% less persistent phase-array storage. It is not a 54.5%
reduction in total process memory. Serial peak RSS(KiB):
M2 115976->113292; M3 113540->113276.
Two-process peaks were114760 and112396KiB per worker; parallel execution
therefore increases total simultaneous process memory.

A fresh review found a pickle regression from storing a ctypes pointer on
the phase object. The pointer was removed. Each process resolves its local
function; a round-trip serialization test now passes. DVODE COMMON remains
protected by process isolation/its existing lock. The new interpolation and
geometry kernels contain no saved state and do not call DVODE.

## Verification

Fresh final integrated suite:73 passed, 8 skipped, 30.56 s, exit 0.
The skipped slow tests are not claimed executed by that pytest run.
The separate complete source campaign above was actually executed.

Native standalone code:
- O0,g,fcheck=all,ffpe-trap=invalid,zero,overflow:31 assertions,exit0.
- O2,ffpe-trap=invalid,zero,overflow:31 assertions,exit0.

Real-equality compiler warnings remain for intentional zero/sign/exact-value
checks; this is not a warning-free build claim.

Additional finite checks:
- existing five-field C19/C20 output 21,360 bytes identical;
- dedicated TOV interpolation against SciPy on n=2,4,333,4097 and 211 queries
  per grid:raw slopes, knot derivatives and dense values identical;
- actual M2/M3 TOV construction with SciPy PCHIP disabled succeeds;
- serialization red/green regression;
- new auditor rejects a corrupted branch status and a nonfinite mass;
- git diff --check succeeds.

One fresh-context read-only review was used; its one actionable finding was
repaired and tested. No review-of-review cycle was run.

## Reproduction

From this worktree:

~~~bash
export PYTHONPATH="$PWD/python_frontend"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

python3 -m pytest -q \
  python_frontend/tests/test_c20_native_tov_pchip.py \
  python_frontend/tests/test_fortran_pchip_geometry_c13_20261121.py \
  python_frontend/tests/test_c15_tov_source_convergence.py \
  python_frontend/tests/test_c17_extended_branch.py \
  python_frontend/tests/test_c17_original_science.py \
  python_frontend/tests/test_c16_source_mode_continuation.py \
  python_frontend/tests/test_cowling_c14_native_error_budget.py

bash scripts/bench_c20_native_tov_pchip.sh 1 /tmp/c20_new_serial
bash scripts/bench_c20_native_tov_pchip.sh 2 /tmp/c20_new_parallel

python3 scripts/audit_c20_tov_migration.py \
  --baseline /path/to/c19_science \
  --current /tmp/c20_new_parallel \
  --output /tmp/c20_new_migration_report.json
~~~

The runner prebuilds the native library once before launching workers, avoiding
the inherited cold-build temporary-name race. This wrapper uses the unchanged
C17 original-source campaign; the scientific schema stays unchanged.

## Remaining scientific and implementation work

Existing Python EOS state evaluation and sample-grid construction remain.
The native pipeline is therefore not yet entirely Fortran behind a thin UI.
No full-GR spectrum, frozen-composition closure, universal root conditioning,
absolute frequency error bound or radial-stability result is asserted here.

The companion COWLING_KINETIC_PRODUCT_MAXWELL_20261123.md derives a physical
same-background kinetic product and ideal slow/rapid interface cancellation.
Its cross-background transport is a diagnostic convention, not a replacement
for C17's acceptance gate. That diagnostic is not implemented in C20.

Next bounded migration: original EOS p/epsilon/n source evaluation and source
grid generation in Fortran, retaining the current source lock and a matched
observable/error budget. The documented raw-candidate roundoff sensitivity
also motivates a separate conditioning/convergence study near M3 onset.
