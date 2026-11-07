# Native Fortran master-function density-gap fill (v1)

## Contract

`compose_eos.fill_density_gap(data, n_fill=60, fit_points=12, gap_phase=None)` selects the largest adjacent density ratio and fills that interval in the native Fortran EOS core. The public Python function only validates the request, serializes the EOSData object, invokes the native executable, and reads the returned EOSData object. It does not call the historic Python numerical routine.

The implementation retains all source rows outside the selected interval. The `n_fill` inserted points are logarithmically spaced in baryon density. For a cold barotrope the master function is the free energy per baryon, with

\[
 x=\ln n_B,\qquad Q_1(x)=\frac{dF}{dx},\qquad
 P=n_B Q_1,\qquad \mu_B=F+Q_1,\qquad \epsilon=n_B F.
\]

Consequently, on the inserted segment,

\[
 dF=Q_1\,d\ln n_B,\qquad
 \epsilon+P=n_B\mu_B,
\]

and the pressure derivative and sound speed are evaluated from

\[
 \frac{dP}{dn_B}=Q_1+\frac{dQ_1}{d\ln n_B},\qquad
 c_s^2=\frac{Q_1+dQ_1/d\ln n_B}{F+Q_1}.
\]

The caller must explicitly pass `gap_phase` when `composition_available` is true and the two source endpoints have different phase IDs. This avoids silently assigning one phase across a boundary. In that case all inserted points receive the supplied phase label.

## Numerical construction

1. A shape-preserving cubic Hermite interpolant (PCHIP) is built for `log(Q1)` against `log(n_B)` using the configured number of nearby source points on either side.
2. The base `Q1` curve is `exp(PCHIP(log(Q1)))`. Its integral over the gap is evaluated by composite Simpson quadrature. A quartic endpoint-flat bump corrects the integral so that the total work equals the source free-energy difference exactly up to quadrature error:

   \[
   \Delta=\frac{F_1-F_0-\int_{x_0}^{x_1}Q_{1,\rm base}(x)\,dx}{x_1-x_0},\quad
   Q_1=Q_{1,\rm base}+30\Delta u^2(1-u)^2,
   \]

   where \(u=(x-x_0)/(x_1-x_0)\). The free energy is obtained by integrating this same `Q1`, not by independently interpolating `F`.
3. The candidate is sampled densely and accepted only when `F`, `Q1`, `mu_B`, `dP/dn_B` and `c_s^2` satisfy the positive, mechanically stable and causal checks.
4. If the corrected PCHIP candidate is not physical, a cubic-Hermite `F(x)` interpolant using `(F0,F1,Q1_0,Q1_1)` is tested from a freshly initialized admissibility state. The fallback evaluates its analytic second derivative, not a duplicate of its first derivative:
   
   \\[
   F''(x)=\\frac{12t-6}{h^2}F_0+\\frac{6t-4}{h}Q_{1,0}
          +\\frac{-12t+6}{h^2}F_1+\\frac{6t-2}{h}Q_{1,1},
   \\qquad \\frac{dP}{dn_B}=F'+F''=Q_1+F''.
   \\]
   
   The fallback is accepted only when `Q1=F'`, `mu_B=F+Q1`, `dP/dn_B`, and `c_s^2=(dP/dn_B)/(F+Q1)` are finite and satisfy the positive, mechanically stable and causal bounds. This explicitly rejects cases where `Q1>0` but `Q1+F''<=0`. If neither construction passes, the operation fails and emits no output table.
5. `Q4`, `Ye`, `Yp`, `Ymu`, `Abar`, `Zbar` and `Yh` use local PCHIP interpolation. `Yn` is set to zero in the inserted cells. The cold-slice conventions `Q2=Q5=0` and `Q7=Q6` are applied to the output, matching the reference Python routine. If both interpolated lepton fraction and nuclear charge are positive, `Abar=Zbar/(Ye+Ymu)` is used to prevent charge overshoot; otherwise `Abar` is interpolated and `Zbar=(Ye+Ymu)Abar`.

The implementation is in `fortran_eos/src/compose_eos_reconstruct.f90`. The CLI action is `fill-gap INPUT_INTERNAL_DIR OUTPUT_INTERNAL_DIR N_FILL FIT_POINTS [GAP_PHASE]`.

## Validation performed in the recovered sandbox

Environment: Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, pytest 9.0.2, GNU Fortran 14.2.0.

- Native Fortran build: passed using `bash fortran_eos/build.sh build/compose_eos_fortran`.
- Frontend and direct-provider tests: **22 passed** using `PYTHONPATH=python_frontend python -m pytest -q python_frontend/tests pipeline/tests`.
- Synthetic constant-`Q1` EOS: inserted `F`, `Q1`, `mu_B`, pressure and energy satisfy the analytic solution; source knots remain unchanged.
- Synthetic `F=930+100 n_B`, `Q1=100 n_B` EOS: reconstructed `F` and `Q1` agree with the analytic solution at relative tolerance `2e-11`.
- Phase-crossing behavior: an unlabelled `0 -> 1` phase crossing is rejected; an explicit phase label is propagated to every inserted row.
- Real source test: `GPPVA_FSU2H` has 300 rows and a maximum adjacent density ratio of 7.4343642318 at source index 216, crossing phase IDs `0 -> 7`. With `n_fill=12, gap_phase=7`, the native operation returns 312 rows. The inserted-segment first-law residual is `2.7643342084e-16` relative.
- Differential comparison against `eos_generator/compose_eos_pkg/compose_eos/reconstruct.py` for the same `GPPVA_FSU2H` input and parameters: maximum relative differences were `8.61e-16` in `n_B`, `8.94e-16` in `Q1`, `2.21e-16` in `Q4`, `1.28e-15` in `Ye`, `1.24e-15` in `Abar`, and below `2e-16` for the other listed fields; `Q3`, `Q6`, `Q7` and `Iphase` matched exactly. This is one representative differential case, not a proof across every EOS table.
- Refusal parity, now covered by an integration test: the `GM_GM1` input with a 191.5518567 adjacent-density ratio is rejected by both the legacy Python routine and the Fortran implementation because neither candidate satisfies the positivity/stability/causality gate. No unchecked EOS is written.
- Hermite fallback regression tests: a case where the integral-corrected PCHIP fails but the endpoint-constrained Hermite master function is physical is accepted; a separate case with `Q1>0` but `dP/dn_B=Q1+F''<0` is rejected. These tests distinguish the true second derivative of `F` from the erroneous duplicate-first-derivative expression.

The full-table cold first-law residual for the GPPVA example remains `1.5379901059e-6` because the original source rows are deliberately retained, including their pre-existing residual. The newly inserted segment itself is consistent to roundoff. This is not a claim that every untouched source knot has been repaired.

## Remaining parity boundary

This increment ports only density-gap filling. It does not port the separate full-table master-function smoothing/candidate-selection pipeline, the six crust-core matching constructions, literature-specific model solvers, genuine Maxwell/Gibbs phase-transition generation, or causal completion. Those remain explicitly unavailable in the Fortran-backed Python API.
