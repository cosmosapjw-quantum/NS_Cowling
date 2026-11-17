# Native Fortran CSS Causal Completion v1

**Status:** implementation and regression tests added to the feature branch; CI verification is the release gate.  
**Branch:** `research/fortran-eos-gap-fill-v1`  
**Target representation:** cold, one-dimensional, strictly monotone CompOSE barotrope  
**Units:** (n_B) in fm\(^{-3}\); (P) and \(\epsilon\) in MeV fm\(^{-3}\); \(\mu_B\) in MeV; \(c_s^2\) dimensionless.  
**Interpretation:** explicit phenomenological high-density continuation, not a microscopic EOS calculation and not a unique repair of the source EOS.

## 1. Scope and relation to the existing Python implementation

The legacy package contains `compose_eos.causality.causal_completion`. It locates the first resolved (c_s^2=1) crossing, matches the pressure and energy density there, checks the baryon chemical potential, and appends a constant-sound-speed (CSS) branch. This document describes a native Fortran implementation of that bounded capability.

The feature is deliberately narrower than the remaining EOS-generator port. It does not port the three-scheme causality profile, full master-function smoothing and candidate selection, crust-core matching, literature RMF/Skyrme/APR solvers, or genuine Maxwell/Gibbs phase-transition construction. In particular, the legacy crust-core token `maxwell` denotes a smooth near-pressure-plateau bridge in that package; neither that method nor this CSS extension should be mistaken for an exact coexistence construction with a density jump.

The CSS parametrization is established in the neutron-star EOS literature. The implementation here uses its constant-(c_s^2) differential closure as a continuation above a numerically identified causal-limit point. This is not identical to every CSS parameterization in the literature, some of which introduce a sharp phase transition and independent transition parameters.

## 2. Thermodynamic derivation

We consider cold, one-component barotropic matter in units where the relation (P=n_Bmu_B-\epsilon) is written without an explicit factor of (c^2). The first law and Euler relation are

\[
  \frac{d\epsilon}{dn_B}=\mu_B,
  \qquad
  \mu_B=\frac{\epsilon+P}{n_B},
  \qquad
  c_s^2=\frac{dP}{d\epsilon}.
\]

Let the matching point be \((n_t,P_t,\epsilon_t)\), and define

\[
  W_t=\epsilon_t+P_t=n_t\mu_t,
  \qquad
  q=\frac{n_B}{n_t},
  \qquad
  s=c_{s,\mathrm{tail}}^2\in(0,1].
\]

The CSS closure is

\[
  \frac{dP}{dn_B}=s\frac{\epsilon+P}{n_B},
  \qquad
  \frac{d\epsilon}{dn_B}=\frac{\epsilon+P}{n_B}.
\]

For (W=\epsilon+P), these equations give

\[
  \frac{dW}{dn_B}=\frac{1+s}{n_B}W,
  \qquad
  W(n_B)=W_t q^{1+s}.
\]

Integration of the pressure equation then yields

\[
  P(n_B)=P_t+\frac{s}{1+s}W_t\left(q^{1+s}-1\right),
\]

and the Euler relation gives

\[
  \mu_B(n_B)=\frac{W(n_B)}{n_B}
             =\mu_t q^s,
  \qquad
  \epsilon(n_B)=\mu_B(n_B)n_B-P(n_B).
\]

Consequently,

\[
  \frac{dP}{d\epsilon}=s
\]

throughout the appended segment. The matching values (P_t), \(\epsilon_t\), and \(\mu_t\) are continuous by construction. The sound speed is continuous at the join only when (s=1), because the matching point is defined by the source-side condition (c_s^2=1). For (0<s<1), the construction intentionally introduces a finite sound-speed kink while retaining continuity of (P), \(\epsilon\), and \(\mu_B\).

The Wolfram symbolic derivative check of the closed-form expressions returns (dP/dn_B=sW_t n_B^s/n_t^{1+s}), (d\epsilon/dn_B=W_t n_B^s/n_t^{1+s}), and their ratio (s).

## 3. Locating the matching point

Set (x=\ln n_B). The source pressure and energy density are represented by shape-preserving cubic Hermite (PCHIP) interpolants of (\ln P(x)) and (\ln\epsilon(x)). On an interpolation segment,

\[
  c_s^2(x)
  =\frac{P(x)}{\epsilon(x)}
    \frac{d\ln P/dx}{d\ln\epsilon/dx}.
\]

The implementation evaluates the derivative at source knots and searches in increasing density for the first bracket with the left value at or below unity and the right value above unity. The interval must satisfy \(\Delta\ln n_B\leq0.25\). A bracket is rejected when its local PCHIP support stencil crosses a change in the source `Iphase` label; the routine does not use a polynomial derivative assembled across such a phase boundary. The crossing is then located by bisection on the PCHIP derivative ratio.

The source must contain at least four points, with positive and strictly increasing density, pressure, and energy density. Non-positive quantities, a non-monotone pressure or energy-density sequence, a mechanically unstable source knot below the first causal-limit crossing, a phase-edge-only crossing, or the absence of a resolved crossing cause an explicit error. A pressure plateau or a density discontinuity is not regularized by this routine.

At the crossing, the PCHIP-interpolated source baryon chemical potential is compared with the Euler value

\[
  \mu_t^{(\mathrm{Euler})}
  =\frac{\epsilon_t+P_t}{n_t}.
\]

The relative difference, normalized by \(\max(|\mu_t|,1\,\mathrm{MeV})\), must not exceed (5\times10^{-3}). This tolerance is inherited from the legacy Python implementation; it is a numerical acceptance threshold, not a universal physical constant.

## 4. Output table semantics

The source rows with (n_B<n_t) are retained in their original density, pressure, and energy-density values. The CSS segment begins with the interpolated matching point and contains `n_points` logarithmically spaced rows through (n_{\max}). The default (n_{\max}) is the maximum density of the source table; an explicit positive (n_{\max}>n_t) may be provided.

The complete output is re-expressed as a cold one-component barotrope:

\[
  Q_1=\frac{P}{n_B},
  \qquad
  Q_3=\frac{\mu_B}{m_n}-1,
  \qquad
  Q_6=Q_7=\frac{\epsilon}{m_n n_B}-1,
\]

with (Q_2=Q_4=Q_5=0). This reconstruction enforces the Euler identity on the output table. It also means that the original chemical-potential and free-energy columns are not copied verbatim, even at retained low-density source knots; the pressure and energy-density values are retained, while the other thermodynamic columns are rebuilt consistently from them.

Composition is intentionally omitted from the entire output, not only from the appended tail. The routine sets `composition_file_present=false`, `composition_available=false`, `nadd=0`, all composition arrays to zero, and all output phase labels to zero. A CSS continuation does not determine the microscopic particle fractions or the nuclear composition. The metadata field `composition_omitted` records whether the source actually carried composition information.

## 5. Interfaces

### Python API

The Fortran-backed frontend exports:

```python
from compose_eos import read_compose, causal_completion, write_compose

source = read_compose("path/to/source_compose")
completed, report = causal_completion(
    source,
    cs2=1.0,
    n_max=None,       # None selects the source maximum density
    n_points=200,
    return_metadata=True,
)
write_compose("path/to/completed_compose", completed)
```

The default return is an `EOSData` object. With `return_metadata=True`, the function returns `(EOSData, report)`, where the report is the native JSON payload containing the match density, pressure, energy density, source and Euler chemical potentials, relative chemical-potential residual, extension sound speed, final density, crossing index, composition-omission flag, and whether the sound speed is continuous at the join.

The API validates the scalar controls in Python and delegates all interpolation, crossing detection, thermodynamic reconstruction, and table generation to the native executable. It does not fall back to the legacy Python implementation.

### Native command line

After compiling with `bash fortran_eos/build.sh build/compose_eos_fortran`, the corresponding internal-CSV command is:

```bash
build/compose_eos_fortran causal-complete \
  SOURCE_INTERNAL_DIR OUTPUT_INTERNAL_DIR 1.0 0 200
```

The fifth argument is (n_{\max}); zero selects the source maximum. The sixth is the number of tail points. The command writes `eos_data.csv` and `eos_data.nml` to the already-existing output directory and emits one JSON status/report line. Non-successful cases emit a JSON error and a nonzero process status.

## 6. Regression strategy

The native tests use an analytic cold EOS

\[
  P(n)=K n^3,
  \qquad
  \epsilon(n)=930n+\frac{K}{2}n^3,
  \qquad K=1000,
\]

for which

\[
  \mu_B(n)=930+\frac{3K}{2}n^2,
  \qquad
  c_s^2(n)=\frac{3Kn^2}{930+\frac{3K}{2}n^2}.
\]

The analytic first causal-limit density is (n_t=\sqrt{930/(3K/2)}\). Tests cover the native API and CLI path, the first-law/Euler identity, source-row pressure and energy preservation, the (s=1) derivative relation, a subluminal (s=0.33) branch, input control validation, rejection when no crossing exists, rejection when the first crossing is tied to a phase boundary, rejection of a source chemical-potential mismatch, and rejection when (n_{\max}\leq n_t).

The native build and all tests must pass in the repository's GitHub Actions workflow. Passing the density-gap tests from an earlier commit is not evidence that this later causal-completion code has passed. Likewise, the existence of this generator does not certify that the resulting table passes every independent EOS consistency test; users should run the native `check_consistency` gate on each generated table and retain its full report.

## 7. Explicit limitations and deferred work

1. The first causal crossing is resolved from source-knot PCHIP derivatives with a maximum accepted knot spacing of (0.25) in (ln n_B). This is not a theorem that excludes arbitrarily narrow between-knot excursions.
2. The source chemical potential is checked only at the match point. A small matching residual does not certify first-law consistency over the entire retained source table.
3. The three-scheme causality profile (PCHIP, natural cubic spline, and centered finite differences) is not yet ported. This routine is not a substitute for that independent diagnostic.
4. The source must be a strictly monotone smooth barotrope. A genuine first-order phase transition, plateau, density jump, or non-monotone segment requires a separate interface-aware treatment.
5. The extension parameter (s), the termination density, and the choice to append a CSS branch are model assumptions. The routine does not infer them from microscopic nuclear theory or observations.
6. Composition is discarded, and the output should not be used as a composition-resolved EOS.
7. This is a bounded port only. Full crust-core matching, all six matching methods, complete master-function smoothing/candidate selection, literature EOS solvers, and the complete phase-transition/TOV treatment remain outside this implementation.

## 8. References

- M. G. Alford, G. F. Burgio, S. Han, G. Taranto, and D. Zappalà, “Constraining and applying a generic high-density equation of state,” *Phys. Rev. D* **92**, 083002 (2015), [doi:10.1103/PhysRevD.92.083002](https://doi.org/10.1103/PhysRevD.92.083002). This paper discusses the constant-speed-of-sound parametrization; its sharp-transition model should not be conflated with the smooth continuation implemented here.
- L. Lindblom, “Causal representations of neutron-star equations of state,” *Phys. Rev. D* **97**, 123019 (2018), [doi:10.1103/PhysRevD.97.123019](https://doi.org/10.1103/PhysRevD.97.123019). This work provides a broader framework for parametrizations that enforce causality and thermodynamic stability.
- M. Fortin et al., “Neutron star radii and crusts: Uncertainties and unified equations of state,” *Phys. Rev. C* **94**, 035804 (2016), [doi:10.1103/PhysRevC.94.035804](https://doi.org/10.1103/PhysRevC.94.035804). This work emphasizes the sensitivity of crust-core matching and the need for causal, thermodynamically consistent matching.
