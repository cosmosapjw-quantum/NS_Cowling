# COWLING-NATIVE-01: smooth polar Cowling DVODE backend

## Source / authority
- Parent commit: `61a78673db57497da1cf991c74cce99a6c2c2c87`.
- Source EOS and central-mass selection: `data/core_support_20261117/primary_v1/REPRESENTATION.json`, `data/stellar_consumers_20261117/tov/CANONICAL_1P4.json`.
- Baseline implementation: `python_frontend/compose_eos/stellar_background.py` and `stellar_cowling.py`.
- Native module: `fortran_eos/src/compose_eos_cowling_dvode.f90`; Python opt-in API: `compose_eos.stellar_cowling.solve_mode(..., backend="fortran_dvode", native_samples=8193)`.
- Main reference remains Python DOP853: default backend is unchanged.

## Physical contract
Metric signature `(-,+,+,+)`, `ds²=-exp(2Phi)c²dt²+exp(2lambda)dr²+r²dOmega²`.
`r` is in km; `omega_g=2 pi f/c_km_s` in inverse km.
For equilibrium-barotropic perturbations with `W=r^(ell+1)w`, `V=r^ell v`, `cs=cs²/c²`:

```text
dw/dr = (omega_g² exp(lambda-2Phi) r v + Phi' w)/cs
        - (ell+1) (ell exp(lambda) v + w)/r
dv/dr = (2 Phi' - ell/r) v - exp(lambda) w/r
```

Integrate in `s=sqrt(h_start-h)`, with `dh/dr=-Phi'` and
`dy/ds=2s/Phi' * dy/dr`. Inner center state `w=1, v=-1/ell`.
Outer state is determined by vanishing Lagrangian surface pressure:
`w=1, v=-Phi'/(omega_g² exp(lambda-2Phi) r)`.
Integrate inner `s=0 -> s_match` and outer `s_surface -> s_match`, then
seek a zero of `det(u_inner,u_outer)/(|u_inner||u_outer|)`.

No frozen composition, solid crust elasticity, dynamical metric perturbations,
stellar rotation, damping times or full-GR quasinormal claims are made.
This is a smooth equilibrium-barotropic Cowling solver.

## Numerical design and limitations
1. The Python exact-contract TOV background is retained as the authoritative
   *input* for this step; background coefficients are sampled once in
   `s=sqrt(h_start-h)`, with shape-preserving cubic PCHIP interpolation.
2. Source `break_h` points are retained with explicit close left/right
   samples to avoid globally smearing finite jumps of the equilibrium
   sound-speed coefficient. This is a numerical approximation, not a
   mathematical interface matching rule.
3. Fortran DVODE Adams `MF=10` and BDF `MF=22` perform the actual Cowling ODE
   evolution. SciPy Brent only brackets the frequency root and forms
   diagnostic profiles from Fortran-returned shots.
4. ODEPACK DVODE uses internal COMMON state. The ctypes adapter serializes
   setup/shots with a process-wide lock. For parallelism, use *processes*,
   not threads sharing DVODE internal state. Native data tables are reloaded
   when changing the active background.
5. The shared library is compiled from the source checkout into a
   content-addressed temporary path. Requires `gfortran` and BLAS
   (`-lblas`). This is an opt-in source-checkout execution path, not yet a
   self-contained binary-wheel installation.
6. An analytic incompressible uniform-density star is used **only** as an
   acausal Kelvin/Cowling limit control; `cs²/c² = +inf` becomes a finite
   very large surrogate for native arithmetic.

**Not implemented:** sharp Maxwell-density-jump slow/rapid conversion, hybrid
star interface matching, multiphase EOS switching or a wholly native TOV
background generator for this adapter. These remain in existing separate
Fortran/Python implementations; this branch does not silently apply a smooth
background kernel to a jump.

## Literature for later interface work
Tonetto & Lugones, *Phys. Rev. D* **101**, 123029 (2020),
https://doi.org/10.1103/PhysRevD.101.123029 : slow-conversion interfaces
require continuity of radial displacement and Lagrangian pressure,
whereas rapid-conversion interfaces obey different junction conditions.
Any subsequent native sharp-Maxwell port must implement and validate both
physical regimes; the current smooth shot has neither.

## Reproduction
```bash
cd /home/cosmosapjw/ns_cowling_wt/cowling_native_dvode_20261118
OPENBLAS_NUM_THREADS=1 PYTHONPATH=python_frontend pytest -q \
  python_frontend/tests/test_native_cowling_dvode_20261118.py \
  python_frontend/tests/test_stellar_cowling.py \
  python_frontend/tests/test_cowling_end_to_end.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=python_frontend \
  python3 scripts/benchmark_native_cowling_dvode_20261118.py
```

Benchmark input, full baseline run cost, 7 native repeats per mode, physical
frequencies, node counts, residuals and code/library SHA are saved to
`data/cowling_native_dvode_20261118/BENCHMARK.json`. Always compare warm
shot speed separately from native setup and compilation cost.


## Verified snapshot on 2026-11-18 (KST)
- Local toolchain: gfortran 13.3.0 on Ubuntu 24.04, Python 3.12.
- Static compilation with `gfortran -c -O0 -g -Wall -Wextra -fcheck=all`
  succeeded. Warnings were limited to unused arguments of the intentionally
  unused DVODE user-Jacobian callback (`MF=10/22`).
- Targeted standalone Cowling + EOS interop tests (including f/p1 real-EOS,
  Adams/BDF parity and hard failure on sharp interfaces): **47 passed**,
  `pytest -q python_frontend/tests/test_cowling_reaudit_contracts.py
  python_frontend/tests/test_css_cowling_adapter.py
  python_frontend/tests/test_stellar_cowling.py
  python_frontend/tests/test_native_cowling_dvode_20261118.py`.
- Supplementary end-to-end subset: **12 passed, 3 skipped**.
- 1.4 Msun exact-contract EOS baseline: radius 11.752078353190377 km,
  source smooth model with 452 interior finite EOS construction breakpoints.
- Python reference (`rtol=1e-9, atol=1e-12`): f=2401.381874932175 Hz,
  p1=6721.957056883109 Hz.
- Native Adams (`native_samples=8193`, 9549 actual knots counting
  one-sided break samples): f=2401.377371448396 Hz, p1=6721.958837881455 Hz.
  Signed errors: f=-0.004503483779 Hz, p1=+0.001780998346 Hz; node counts
  (W,V)=(0,0) and (1,1). The matching residual magnitudes are
  1.22e-7 and 5.43e-8.
- Independent native BDF verification: f=2401.377332881328 Hz,
  p1=6721.962385534773 Hz. Adams/BDF disagreement <0.004 Hz for both.
- Final same-host benchmark after rebuild: single-run Python DOP853
  f=19.680 s, p1=28.444 s; median of seven warm native runs
  f=0.02639875 s, p1=0.03272017 s. Ratios approx **745.5x** and
  **869.3x**, respectively. These are warm shooting comparisons,
  not total source loading, native compilation, TOV background, or
  many-background production throughput. Native input sampling alone took
  approximately 0.5 s on this host.
- Native 16385-sample cross-check: f=2401.3812598631125 Hz,
  p1=6721.959660883174 Hz, showing stability better than 0.01 Hz
  relative to the 8193-sample results for this star.
- All numeric benchmark inputs and executable source hashes are in
  `data/cowling_native_dvode_20261118/BENCHMARK.json`.

## Remaining scope
**Status:** smooth Cowling shooting backend IMPLEMENTATION_VERIFIED and
bounded numerical parity checked for one realistic reference star plus
one Kelvin weak-gravity limit. Not yet production-grade native whole-star
pipeline. `COWLING-NATIVE-02` should replace the Python TOV/background sampling
with native Fortran and perform mass/EOS-grid reproducibility checks.
`COWLING-MAXWELL-NATIVE-03` should port slow and rapid sharp-interface
jump conditions with independent junction residual controls; do not
accidentally identify their eigenproblems.
