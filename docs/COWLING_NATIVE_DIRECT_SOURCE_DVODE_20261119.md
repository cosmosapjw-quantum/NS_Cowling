# NS_Cowling: physical source–evaluating DVODE and sound-speed knot restarts

2026-11-19 KST · Stage C5 (implementation), C6 (regression), C8 (accuracy).
Baseline: \`daf5de300b06009f8654f795aa99fe40d512666e\`.
Branch: \`research/cowling-adaptive-coefficients-20261119\`.
Native host subagents: \`NOT_RUN\` (no delegate/spawn facility exposed);
the parent executed a serial bounded scientific coding task.

## Scientific question and actual defect

Previously the source-locked SPG(M3)+PK1 CSS 0.7 solar-mass slow
\`p1\` mode was **not admitted** by the original Fortran DVODE method.
It was admitted by a Python DOP853 oracle, with f=3736.208895 Hz.
On the 1.4 solar-mass W1/V2 mode, the old native answer changed by
~0.609 Hz between 4097 and 8193 Cowling coefficient samples.
The old native solver did converge for the *interpolated ODE*, but
the source coefficient interpolant differed from the original
one-sided barotropic source.

The new physical problem is not solved by simply reducing DVODE \`rtol\`,
raising the PCHIP coefficient grid uniformly or declaring a smaller
shooting residual to be a smaller eigenfrequency error.

### Measured source interpolation pathology

The code uses \`s=sqrt(h_start-h)\` and the TOV background state
coefficients \`(r, Phi, lambda, Phi_prime, cs2/c2)\`. Original generated
piecewise EOS internal knots have *discontinuous derivative* \`cs2\`
even within one physical phase. For the existing M3 CSS 0.7 solar-mass
source and PCHIP-4097 setting, evaluating every intermediate grid midpoint
against its original **one-sided EOS** gives:

| Physical phase | Native PCHIP coefficient grid | Original EOS internal breakpoints | Selected strong sound-speed jumps (>1%) | Max pointwise normalized cs2 error |
|:--|--:|--:|--:|--:|
| High-density | 6601 | 1065 | 1 | 0.32584 |
| Lower-density | 2710 | 332 | 39 | 189.4153 |

Near one of the low-density knots, \`cs2_exact≈0.00263275\`
but the unrefined interpolation produces \`cs2≈0.501315\`
within a **~10^-10-wide s gap**. This enormous *local coefficient
comparison* does **not** imply a 189-fold eigenfrequency error;
it proves that source derivative jumps cannot be certified
by assuming globally smooth PCHIP coefficients.
Source evidence: \`data/cowling_direct_source_callback_20261119/KNOT_AUDIT.json\`.

## Code changes

Two source files carry the physical numerical change:

- \`fortran_eos/src/compose_eos_cowling_dvode.f90\`
- \`python_frontend/compose_eos/stellar_sharp_native.py\`

An existing ctypes signature file additionally changes:
\`python_frontend/compose_eos/stellar_cowling_native.py\`.

All code works in the original (-,+,+,+) convention, with \`c\`
maintained in \`omega_geom=2*pi*f/C_KM_S\`.
For the native independent variable
\[
s=\sqrt{h_0-h},\qquad
dh/ds=-2s,\qquad dr/dh=-1/\Phi'(r),
\]
so the source-evaluated original Cowling ODE obeys
\[
\frac{d\mathbf y}{ds}
= \frac{2s}{\Phi'(r)} F_r(\mathbf y;\,r,\Phi,\Lambda,\Phi',c_s^2).
\]
The source lookup uses the **same original EOS and TOV** whose
phasewise interface physics was already audited, never a fitted
or smoothed substitute.

1. **Opt-in direct source backend:** Python registers a process-local
   \`ctypes.CFUNCTYPE\` callback into the Fortran DVODE module.
   At each actual RHS evaluation the callback returns the original
   one-sided \`state_at_enthalpy\` coefficients. The deliberately sparse
   33-point phase table is only a native endpoint sentinel, **not**
   the coefficient interpolant. \`samples=129\` in this mode does
   NOT mean 129 physical RHS samples: thousands of actual source calls
   occur during shooting. The callback is serialized by the existing
   native global lock, strongly retained until the F90 call returns,
   and unregistered in \`finally\`. Exceptions become Fortran
   status failures and preserve the original error type/message.

2. **Restart on genuine strong source EOS derivative jumps:** Evaluate
   the original \`c_s^2\` on both sides of each *internal* EOS
   \`break_h\` in a fixed narrow physical neighborhood
   \`dh=max(2e-11, 1e-9*max(h,1e-5))\`.
   Only \`|cs2_hi-cs2_lo|/max(cs2_hi,cs2_lo)>=0.01\` receives
   a true fresh DVODE restart. The M3 low-mass example has 1 high
   phase plus 39 lower phase strong knots (40 total), rather than
   restarting at all 1397 original phase-local EOS breakpoints.
   The F90 \`cowling_native_segment_restarts\` transports state
   exactly through each selected source knot, restarts DVODE history
   at the event, and preserves profile output requested by the caller.

3. **Correct one-sided limit at restart:** At the beginning or end of
   an integrator subinterval, the F90 source callback evaluates the
   derivative from the **current** side, at an interior offset
   \`min(1e-9,0.01*Delta_s)\` in s. This avoids choosing an arbitrary
   derivative value at an EOS \`C^0\` knot. It does not smear a
   physical Maxwell density jump: the latter remains a separate
   phase with the originally implemented slow/rapid transfer.

The original default remains
\`NativeSharpCowlingDVODE(...,coefficient_backend="pchip")\`;
no old API arguments are altered. The new documented alternative is:

\`\`\`python
native = NativeSharpCowlingDVODE(
    background, CowlingSettings(surface_h=3e-8, rtol=1e-9,
                                atol=1e-12, max_step_h=0.002),
    samples=129,
    coefficient_backend="direct_source",
    source_restart_relative_jump=0.01,
)
mode = native.solve_mode((3710,3760),conversion="slow",profile_points=350)
\`\`\`

An explicit \`source_restart_relative_jump=None\` disables internal
source derivative restarts but still uses direct source callbacks.
This is intended as a falsifiable diagnostic, not a recommended
high-accuracy production setting.

## Actual physical numerical outcomes

All \`Delta f\` values compare the recomputed mode against the
**separately executed** source-EOS Python DOP853 baseline with its own
TOV integration at the same original EOS and center density.
They are observations, **not** a rigorous a posteriori error bound.

| Source-locked star | Mode, interface condition | Old native PCHIP-8193 | New direct-source DVODE after knot restart | Direct-source minus DOP853 |
|:--|:--|:--|--:|--:|
| SPG_M3+PK1 CSS, 0.7 Msun | slow, strict W1/V1 | Rejected, residual >1e-5 | **3736.210476 Hz** | +0.001581 Hz |
| SPG_M3+PK1 CSS, 1.4 Msun | slow W1/V2 conditional | 5329.147252 Hz | **5329.202223 Hz** | −0.000451 Hz |
| SPG_M2+PK1 Maxwell, 1.4 Msun | slow, strict W1/V1 | 5147.989067 Hz | **5148.000868 Hz** | −0.000162 Hz |
| SPG_M2+PK1 Maxwell, 1.4 Msun | rapid, strict W1/V1 | 5143.704239 Hz | **5143.719701 Hz** | −0.0000206 Hz |

The old native PCHIP-16385 original coefficients in this benchmark
are still insufficient to reproduce the W1/V2 M3 result
within 0.03 Hz; the old default 8193 and 16385 PCHIP
calculations **remain numerically unchanged** before and after
this implementation. A lower residual on its own would not have
detected the source error.

The rapid mode has source-DOP853 f=5143.719722 Hz,
interface traction residual ~4.37e-16 in the native calculation
and relative join ~3.5e-7. This is a new actual rapid test and
does not alter either conversion condition.
Slow mode M3 at 1.4 Msun remains **W1/V2**, never renamed strict p1.

### Sensitivity to DVODE tolerance after one-sided restarts

Same source M3 CSS 0.7 Msun, same native TOV and slow conversion,
with the following mode integrator changes only:

| rtol | max_step_h | f (Hz) | relative join | admission |
|:--|--:|--:|--:|:--|
| 1e-8 | 0.0020 | 3736.210052 | 3.28e-6 | PASS |
| 1e-9 | 0.0010 | 3736.210460 | 5.91e-8 | PASS |
| 1e-10 | 0.0005 | 3736.210542 | 2.31e-6 | PASS |

Observed frequency spread: ~0.000490 Hz.
**Prior to restart**, only the middle of these three
independently executed configurations passed the fixed
native physical admission threshold. The exact before/after
receipts are separate: \`TOLERANCE_BEFORE_RESTART.json\` and
\`TOLERANCE.json\`. Tightening \`rtol\` alone does not
guarantee monotone decrease of the relative matching residual.

### Measured wall times (non-general performance)

From one workstation/process, with thread counts locked to 1,
each **same-EOS, same-mode** wall time is reported by the script,
including its own native phase setup/root shooting. For M2
Maxwell 1.4 Msun rapid p1:

- original PCHIP-8193: **1.32 s**, \`|delta f|≈0.01548 Hz\`;
- new direct-source F90: **5.92 s**, \`|delta f|≈2.06e-5 Hz\`;
- independent Python DOP853: **32.01 s**, direct reference.

Direct-source F90 is slower than PCHIP at lower resolution, but in
this one scenario it is around 5.4x faster than the DOP853 source
reference. It performs ~84k exact Python state callbacks over many
frequency shots, not just 129 coefficient evaluations.
No matched-accuracy wall-time comparison with a fully converged
PCHIP reference is available across the entire EOS family.

## Claim admission, controls and remaining physics

- \`implementation-verified\`: new explicit C callback contract,
  safe function pointer lifetime/reset on failure,
  source SHA-locked derived EOS; existing native executable
  API and original default unchanged; one-sided phase restarts.
- \`numerically-checked\`: physical 0.7/1.4 Msun slow and
  1.4 Msun rapid bounded source-vs-PCHIP-vs-DOP853 campaign.
  The 0.01 relative jump threshold is a **local pragmatic
  numerical selection rule**, not a microphysical transition threshold.
- \`unresolved\`: formal global frequency error bounds and
  general optimal knot selection; globally continued avoided-crossing
  mode family; all 96 stars and 192 spectra; original author SPG/PK1
  microscopic crust-core unification; elastic crust, superfluid
  composition-frozen perturbations and full metric GR modes.
- \`not changed\`: physical slow \`[xi^r]=0\` and rapid
  \`[xi^r-Delta P/P_0']=0\`, together with \`[Delta P]=0\`
  for both limits. Reference: Tonetto & Lugones
  (2020), Phys. Rev. D 101, 123029
  https://arxiv.org/abs/2003.01259 .
  F90 interface transfer and TOV source construction remain fixed.
- \`no publication gate\`: a small bounded verification sample
  cannot be silently promoted to a general-purpose convergence theorem.
  A hybrid result with W1/V2 is a conditional numeric family,
  not a strict p1 label.

**Next bounded optimization opportunity:** avoid Python state callbacks
at every single F90 RHS call by generating phasewise source-aware
quadratures/cached coefficient basis aligned with the documented
strong derivative jumps, then benchmark at *matched* frequency
accuracy and retain the direct-source route as independent control.

## Reproduce

The fresh scripts, immutable exact-density source-locked input,
accepted/rejected statuses and all benchmark receipts are in
\`data/cowling_direct_source_callback_20261119\`.

\`\`\`bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=python_frontend
python3 scripts/audit_sharp_sound_speed_knots_20261119.py
python3 scripts/benchmark_sharp_direct_source_callback_20261119.py
python3 scripts/check_direct_source_tolerance_20261119.py
python3 scripts/verify_direct_source_rapid_20261119.py
pytest -q python_frontend/tests/test_sharp_cowling_exact_source_callback_20261119.py \
  python_frontend/tests/test_sharp_cowling_direct_source_evidence_20261119.py --run-slow
\`\`\`

The default Fortran library compiled by
\`compose_eos.stellar_cowling_native._native_library\`
is keyed by exact Fortran source SHA-256; the original
raw CompOSE GPPVA/secondary BBP/nucleardatapy binaries
are local untracked inputs, not republished.
