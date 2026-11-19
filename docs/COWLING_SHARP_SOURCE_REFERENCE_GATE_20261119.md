# COWLING-SHARP-SOURCE-REFERENCE-GATE-06

Date 2026-11-19 KST. Authoritative parent:
1e7823ecd62a2914d40c0948d3811410b1ea16da.
Worktree branch: research/cowling-sharp-certified-modes-20261119.
Native delegated subagent capability: NOT_RUN (no host tool exposed);
scientific coding performed by the parent in serial fallback.
Requested coding stages: C5 implementation, C6 implementation verification,
C8 numerical accuracy. C2 source/code correspondence was inspected
as necessary; no independent claim that the full ten-stage cycle ran.

## The specific physical/numerical defect

The exact same SPG(M3)+PK1 CSS 1.4 Msun sharp barotrope has
a well-identified W=1, V=2 slow eigenfrequency near 5329.20 Hz.
At a native PCHIP background coefficient grid of 4,097 samples,
the old Fortran DVODE solver finds approximately 5328.5383 Hz and
a **very small** matching residual around 5.4e-7.
At 8,193 samples it finds approximately 5329.1474 Hz,
also with a small residual. This ≈0.6091 Hz shift is
**not a root-solving tolerance error**: it reflects different
interpolated coefficient functions.

At M3-CSS 0.7 Msun the old native p1 was rejected by the physically
fixed matching residual threshold 1e-5, including when increasing
coefficients from 4097 to 8193 samples, while independent phasewise
DOP853 source-EOS evaluation returns the same regular eigenmode.

The objective is to make that accuracy and admission boundary
machine-readable, not to loosen physical tolerances or claim the
numerical method now converges faster. No change to the original
sharp EOS, conversion conditions, Fortran DVODE code,
or pre-existing single-mode API is made.

## Mathematical source, conventions and admission

Underlying field fluid approximation is the established
relativistic Cowling equation in a static (-,+,+,+) background,
with c carried in omega=2*pi*f/c. Bulk closure is
Gamma1 = Gamma_eq (one barotropic degree of freedom),
and there are no metric, elastic or dissipative perturbations.

At a sharp phase boundary the physics conditions are
[Delta P]=0 for both conversion limits; slow additionally
[xi^r]=0, rapid additionally [xi^r-Delta P/P0']=0.
These are the conditions physically derived for
radial/nonradial hybrid-star interfaces in:

- Tonetto & Lugones, Physical Review D 101, 123029 (2020),
  https://arxiv.org/abs/2003.01259
- Pereira, Flores & Lugones, Astrophysical Journal 860, 12 (2018).

The existing native interface transfer is preserved;
this stage only changes how independently obtained
eigenfrequencies are admitted and represented.

To call a native result "crosschecked" we require:

1. Two consecutive actual attempted native coefficient
   resolutions both produce root-admitted eigenfunctions,
   with normalized join residual <=1e-5;
2. The same W and V phasewise node pair, frequency step
   |Delta f|<=0.05 Hz for the bounded campaign;
3. Phasewise sign-invariant overlap of (W,r*V) >=0.98,
   comparing one physical phase at a time;
4. At the *same unchanged physical EOS* and center density,
   direct EOS DOP853 integration yields independent
   f within 0.05 Hz and matching residual <=1e-6;
5. A strict target label requires (0,0) for f or (1,1)
   for p1. The (1,2) W1/V2 branch is explicitly
   CONDITIONAL_NUMERIC_FAMILY_NOT_STRICT_P1.

These criteria are **not** formal error bounds, numerical
proof of global eigenvalue continuation, or microphysical
phase-transformation laws. They define scoped numerical
agreement. The original 1e-5 native admission is never weakened.

## Result: three real, source-locked stars/modes

The original derivation input is the immutable previously
published stage 05 7-ZIP:
data/eos_missing_data_recalc_20261119/artifacts/
EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip.
Each one-phase barotropic 3-layer EOS contract is extracted
from the stored 05_SPG5_PK1_10_WHOLE_STARS.zip and checked
against its internal SHA-256 manifest; no raw author source is
relabelled as reconstructed original SPG or microscopic PK1.

| Model | M [Msun] | Conversion | f [Hz] adopted | Native coefficient data | Admission / node status |
|:--|--:|:--|--:|:--|:--|
| SPG_M3+PK1 CSS | 0.7 | slow | 3736.208895158 | 4097 & 8193 roots rejected | Direct DOP853 only, strict (1,1) |
| SPG_M3+PK1 CSS | 1.4 | slow | 5329.202673892 | 4097: 5328.538310; 8193: 5329.147443 Hz | Direct DOP853 only; W1/V2 NOT strict p1 |
| SPG_M2+PK1 Maxwell | 1.4 | slow | 5148.004699554 | 8193: 5147.989212; 16385: 5148.004700 Hz | Native against direct DOP853, strict (1,1) |

The direct DOP853 reference yields an independent
M2 Maxwell frequency of approximately 5148.001029 Hz;
the native/local reference difference is 0.0036702 Hz
and the native-grid frequency difference is 0.0154874 Hz.
The two native mode shapes have absolute overlap
0.999999998894.

For the M3 CSS 0.7 Msun reference, DOP853 match residual is
3.85e-10. Native resolution 4097/8193 instead gives
approx 6.92e-5/2.26e-5 join residual, both above
the old 1e-5 gate.

For M3 CSS 1.4 Msun the DOP853 one-sided node pair is (1,2)
and the reference join residual is 7.48e-10.
The native f shift from 4097 to 8193 is 0.609132 Hz,
far above the 0.05-Hz campaign admission.
**Do not** label the W1/V2 branch unconditionally p1
on the basis of matching frequencies alone.

The three cases ran in roughly 47/37/39 wall seconds
respectively. A direct source DOP853 reference is
computationally costly and should be scheduled for
selected anchors, not misrepresented as an optimized
production mass-sweep replacement.

Machine-readable source + output:
data/cowling_certified_modes_20261119/RESULT.json.
Every case contains source-contract sha256, native
attempts, physical rejection reasons, reference mass,
frequency, node category, join residual and numerical gates.

## Test and falsifiability

A new real-star test with --run-slow exercises the problematic
M3 CSS 0.7 Msun case from source-locked 05 ZIP.
A separate M3 CSS 1.4 test asserts that a (W,V)=(1,2)
frequency is never promoted to strict p1. Pure tests verify
the no-false-PASS behavior on rejected intermediate grids,
misaligned phase interface signs and changed node counts.

An actual red-green unit experiment on native coefficient
sampling revealed a defect in the initial gate: two accepted
endpoint sample grids could accidentally skip a rejected
intermediate resolution. The fix now requires **two consecutive
attempted grids**, with the last attempted grid also accepted.
This preserves the empirical failure instead of cherry-picking
successful frequencies. Original slow/rapid matching code
remains unchanged.

## Claim gate

- Source/data: previously published source-derived contracts,
  SHA verified; original CompOSE GPPVA and source SPG EDF
  identity not upgraded here.
- Numerical: bounded direct DOP853-vs-native agreement, with
  an explicit native failure and conditional W1/V2 branch.
- Implementation: accepted only after fresh run of the
  focused real tests and prior TOV/Cowling regression suite.
- Scientific: **no universal 1e-7 Hz precision claim**, no
  96-star/192-spectrum full-sweep completion, no radial
  stability theorem and no original-field/generator equality.
- Next coding stage (only if selected): an adaptively refined
  native background and Cowling coefficient representation that
  actually reduces coefficient interpolation errors at matched
  frequency accuracy without always using Python DOP853.
