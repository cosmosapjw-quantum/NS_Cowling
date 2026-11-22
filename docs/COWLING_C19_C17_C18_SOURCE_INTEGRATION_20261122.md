# NS_Cowling C19: source-locked C17/C18 native Fortran integration and scientific replay

Date: 2026-11-22 KST. **This is one bounded implementation/verification
phase, not a new physical mode or global scientific release.**

## Source identity, claims, and immutable boundaries

- C17 source-locked branch: `c4a777013aed52208e3736ead8de76611f7bf4cc`,
  `research/cowling-c17-extended-mode-branch-20261122`.
- C18 isolated improvement: `323bc2b369e6007dd61e4f851fbde39620ce3790`.
  Cherry-picked onto C17 in the C19 worktree as `d4c964d02f588e620da10885d80dd1ed4538007a`.
- New isolation: `research/cowling-c19-source-locked-overlap-integrated-20261122`
  under `/home/cosmosapjw/ns_cowling_wt/cowling_c19_c17_c18_integration_20261122`.
  No update to C17's worktree or default branch.
- Original EOS archive:
  `data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip`,
  byte identity SHA-256
  `299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04`.
  Its M2, M3 original-source contract fingerprints stay unchanged.
- Metric signature (-,+,+,+); original EOS units MeV/fm^3, geometric
  TOV km and km^-2; f in Hz and omega_geom=2*pi*f/(299792.458 km/s).
- Exactly one sharp Maxwell interface; M2 rapid conversion and M3 slow
  conversion; original source EOS derivatives/c_s^2 and cold
  equilibrium Gamma_1=Gamma_eq. TOV sample grid 2049, Cowling 4097,
  original 1e-5 normalized determinant and 1e-5 full-state join gates
  unchanged. No full-GR metric perturbations or frozen composition.

## Observed, independently reproduced C18-to-C17 integration regression

On the unmodified C17 source, the C18 checked Fortran test
`cowling_c18_overlap_stability.f90` failed RED with exit 14, because
a tiny but finite overlap underflowed. The cherry-picked C18 fix
passed its original finite-amplitude synthetic contract.

However the FIRST actual original-source 40-star C17 replay after
cherry-picking C18 failed source-locked branch continuation:
the frequency scans remained BYTE-identical to original C17, but
40-star f/p1 branch acceptance collapsed to an anchor and typed
SHAPE_OVERLAP_REJECTED at +/-1% density. This was NOT a change to
the stellar shooting RHS or to the EOS. The cause was the new
C18 requirement that all radii strictly increase within a phase.

The physical C17 mode profile has 240 samples in the high phase,
120 in the low phase. The two high-phase ODE shooting legs both
include the *same interior matching radius* (high-phase indices
119/120 in zero-based local indexing). Two independently integrated
W values differ slightly there even when their independently gated
physical root and Maxwell join are admitted. Concrete observations:

- M2 @ density factor 0.99, f mode: r=6.151518467347029 km twice,
  W=0.03394988348299955 / 0.033949885963253275.
- M3 @ factor 1.00, f mode: r=5.132184501348457 km twice,
  W=0.015360885979838146 / 0.015360882073162616.

This radius duplication is a legitimate, numerically matched
*shooting interface inside a single phase*, not the true Maxwell
energy-density jump. Treating it as invalid was a numerical-mesh
contract error; it should not be called a new physical instability.

The RED native Fortran C19 test demonstrates this exact case:
`cowling_c19_shooting_join_mesh.f90`, exit 10 and message
`C19_FAIL_C17_SHOOTING_JOIN_REJECTED 1`.

## Narrow Fortran repair and mathematical scope

The only production edit beyond the C18 import is
`fortran_eos/src/compose_eos_cowling_mode_overlap.f90`.

For each phase P, permit **at most one** repeated areal-radius
sample r[j]=r[j+1] when

`|W[j]/max_P|W| - W[j+1]/max_P|W|| <= 1e-5`.

The duplicate amplitude guard uses the same numerical magnitude
as the original 1e-5 join requirement but a *different* phasewise
normalization. It is a separate representation check, and does not
replace the independent root/interface/full-state join gates.
Nonfinite input, decreasing radii, a second internal duplicate or
a large W discontinuity still fails closed. The true Maxwell phase
boundary may repeat radius without forcing W continuity, because
the physical slow/rapid transport laws govern that boundary.

Existing C18 phase-local and global normalized overlaps remain
scale-stable at W magnitudes 1e-200 and 1e200. Interpolation remains
strictly within each physical phase. The 65-point equal-phase-weight
trapezoid over dimensionless areal radius u remains a descriptive
W-only shape correlation, NOT a conserved/covariant pulsation inner
product. No new physical acceptance gate was introduced.

## Actual 40-star reproduction with original physics

Command, 2 independent ODEPACK-safe OS processes:

```bash
bash scripts/bench_c17_extended_parallel.sh 2 /tmp/ns_c19_20261122/source_replay_fixed_2proc
```

This used the unmodified C17 source-exact 111-point frequency grid
and all 20 central densities for each M2/M3 original source.
One- and two-process calculations were independently repeated:

| Process workers | Actual wall seconds | Worker peak RSS (KiB) |
|---|---:|---:|
| 1 | 66.168498 | 115392 / 112876 |
| 2 | 34.770320 | 115032 / 113708 |

Speedup by independent OS process isolation for these two stars:
**1.903x**, not a single-star GPU/OpenMP or DVODE thread-speed claim.
A final *exact-source-byte* independent two-process replay after
the last comments-only source edit used source SHA256
`d7a013e522ca79ca923e4dea34ede6d3714e2a292d0c767d9d32937640d384bd`,
measured 31.826131 s and generated the exact same scientific hash.
Timing differences between 2-worker repeats are runtime variation,
not a demonstrated further optimization or fixed throughput ratio.
The *full scientific output* for 1 and 2 processes was byte identical:
`SCIENTIFIC.sha256` SHA-256
`dd033582b765c31828165d5dd9c83670de256f1610d539f765e87a1409e90254`.

A separate frontend-only, exact-identity/typed-difference tool
`scripts/audit_c19_original_source_replay.py`
checked original C17 versus C19 complete M2 and M3 JSON and raw
original-source scans, without implementing any physics in Python.
It reports `PASS_SCOPED_C17_C19_REPLAY`:

- Both M2 and M3 compressed 2220-sample original Fortran scan
  payloads BYTE-identical to canonical C17 (4440 scan points total).
- All EOS source fields, stellar backgrounds, Maxwell junction/root
  outputs, frequency values and physical branch statuses exact.
- M2: 20/20 f, 16/20 p1 branch admissions unchanged.
- M3: 16/20 f, 7/20 p1 branch admissions unchanged.
- 86 serialized W-overlap floating-point values differ only by normal
  arithmetic order, maximum absolute difference 5.551115123125783e-16.
  Threshold 0.95 and all C17 rejection statuses unchanged, including
  M3 @ factor 1.05 p1 high-phase overlap below threshold.
- Frontend report:
  `data/cowling_c19_integrated_20261122/SEMANTIC_C17_C19_REPLAY.json`.
- Auditor negative-control: a locally corrupted C19 M2 f-branch status
  (`INCORRECT_ADMISSION_INJECTED`) forced exit 1 and
  `FAIL_REPLAY_MISMATCH`; the comparator is demonstrably red-capable.
  This control was run outside the source tree and does not alter the
  canonical 40-star output.

The Fortran kernel-only benchmark `scripts/bench_c18_fortran_overlap.sh`,
against the exact C16 pre-optimization implementation at equal -O2,
n=360 synthetic 2-phase profiles and 20k comparisons with 5
alternating repeats, gave median:

| Workers | Original C16 kernel seconds | C19 kernel seconds | Kernel-only speedup |
|---|---:|---:|---:|
| 1 | 0.1397669 | 0.0820026 | 1.7044x |
| 4 OpenMP | 0.0375916 | 0.0220538 | 1.7045x |

This is an **isolated reentrant overlap kernel**. The larger real-star
computation uses 1/2 isolated processes for non-reentrant DVODE,
with `OMP_NUM_THREADS=1`; it was NOT accelerated by four
OpenMP threads. A full-pipeline gain over C17 itself is NOT claimed.

## Fresh validation and known limitations

Native Fortran `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`
and `-O2 -ffpe-trap=invalid,zero,overflow`: both modes passed
C18 finite-amplitude/mesh regression, C19 actual shooting-join
continuity/large-mismatch/multiple-duplicate regression and
C17 source-onset/branch-selector test (6 program runs).
Legacy test's pre-existing continued-string warning is
not a native numerical failure.

Focused actual source C14-C17 integrated Python frontend/native
regression: **40 passed, 6 skipped, exit 0**.
The full historic 166-test suite was not rerun in C19.
No proof of a self-adjoint/covariant mode inner product,
absolute eigenfrequency error bound, global avoided-crossing
absence or near-onset radial stability follows from this replay.
M2 onset remains outside sampled low-density coverage.

Literature/physical-scope notes:
- Lindblom (1998): https://arxiv.org/abs/gr-qc/9802072
  establishes M-R curve phase transition geometry, not this
  W-only inner product.
- Tonetto and Lugones (2020): https://arxiv.org/abs/2003.01259
  derive slow/rapid phase conversion conditions and discontinuity g
  modes, but do not license this C19 correlation as a conserved product.
- Dong, Kuan and Yunes (2026), arXiv:2609.33002, Phys Rev D
  **accepted 19 November 2026**:
  https://journals.aps.org/prd/accepted/10.1103/r5vk-95v6
  examines f/p1 on full-GR, SLy4+CSS backgrounds and is a relevant
  *future* comparator, not direct EOS-matched validation of this
  equilibrium Cowling code.

**Next bounded physics action:** near-onset matched-EOS high phase
regularity, physical mode metric/inner product and a physically
motivated branch identity criterion, with the current source lock
and explicit TOV/frequency numerical uncertainty maintained.
Do not mutate existing C17/C18 canonical science records;
retain this integration as a new branch and backup.
