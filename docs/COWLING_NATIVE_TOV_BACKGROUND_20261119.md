# COWLING-NATIVE-02: Fortran DVODE TOV background + Cowling modes

## Scope and source identity

Parent commit: `e0f481178b327ff344b15152679a5b192c16ad90`, which
completed the earlier smooth Cowling-only native shooting stage.
Current changes are isolated to branch
`research/cowling-native-tov-background-20261119`, not the dirty main
NS_Cowling worktree. No automatic merge is performed.

Goal: replace the SciPy TOV *ODE integrator* with Fortran DVODE and deliver
its radial interior profile to the already native Cowling shooting kernel.
Retain the original exact-contract EOS as a read-only Python source
provider for the one-time sampled enthalpy table; do **not** claim a
completely Python-free EOS implementation.

Sources of scientific reference:
- `python_frontend/compose_eos/stellar_background.py`, exact piecewise EOS
  and SciPy DOP853 TOV;
- `python_frontend/compose_eos/stellar_cowling.py`, SciPy DOP853 Cowling;
- `data/core_support_20261117/primary_v1/REPRESENTATION.json`;
- `data/stellar_extension_20261117/eos/{medium,fine}/REPRESENTATION.json`.
- Independent control: analytical Gamma=2 relativistic cold polytrope,
  `P=K n²`, `epsilon=m_b n+K n²` with K=400 and m_b=939 MeV,
  defined in `python_frontend/tests/test_native_polytrope_cowling_20261119.py`.

The supplied medium and fine representations share substantial overlapping
underlying physics with primary_v1; treat them as numerical/
representation variants, **not independent nuclear-matter models**.

## Continuous equations and convention

Metric signature (-,+,+,+):
`ds²=-exp(2Phi)c²dt²+exp(2lambda)dr²+r²dOmega²`,
`exp(-2lambda)=1-2m/r`. The enclosed geometric mass `m=GM/c²` and radius
`r` are in km; pressure `P` and energy density `epsilon` in km^-2
after source conversion `G P_SI/c^4 × 10^6`.
The original enthalpy is `h=∫dP/(epsilon+P)`.

In enthalpy h, the TOV equations are:
```text
dr/dh = -r (r-2m)/(m+4π r³P)
dm/dh = 4π r² epsilon (dr/dh)
dPhi/dr=(m+4π r³P)/(r(r-2m)).
```
For `s = sqrt(h_start-h)`, with `h=h_start-s²`:
```text
dr/ds = 2s r (r-2m)/(m+4π r³P)
dm/ds = 4π r² epsilon (dr/ds).
```
`s=0` starts from the regular small-center-radius expansion at
`r0=0.001 km`, `m0=4π epsilon_c r0³/3`,
`h_start=h_c-2π(epsilon_c/3+P_c)r0²`. The approximation order is
the same as the Python reference, not a higher-order center expansion.

For a cutoff `h_surface=1e-4`,
`Phi(h)=0.5 log(1-2M/R)-h+h_surface`,
`lambda=-0.5 log(1-2m/r)`.
Use the identical Cowling `ell=2` equation and matched center/surface
conditions recorded in `docs/COWLING_NATIVE_DVODE_20261118.md`.

## Code surfaces and execution

- `fortran_eos/src/compose_eos_cowling_tov_dvode.f90`: Netlib DVODE
  integration of r,m at requested s points, using in-module PCHIP EOS
  coefficients and physical domain validation. `ITASK=4` and
  `RWORK(1)=s_surface` prevent RHS evaluation outside EOS support.
  `MF=10` Adams is default, `MF=22` BDF supported.
- `python_frontend/compose_eos/stellar_cowling_tov_native.py`:
  immutable-source EOS state sampling, center initial values, native
  arrays, metric reconstruction, vectorized Cowling background coefficient
  accessor, and provenance. Python performs no TOV integration.
- `python_frontend/compose_eos/stellar_cowling_native.py`: native
  library now includes both TOV and Cowling modules. The native TOV object
  feeds vectors without Python `solve_ivp` calls. Preserves Python-only
  Cowling API default.
- `python_frontend/compose_eos/stellar_background.py`: optional
  `backend="fortran_dvode"` and `native_samples=8193` on
  `solve_background`. Default `backend="python"` is unchanged.

A shared process-wide lock protects the non-reentrant DVODE COMMON state.
For parallel operation, run separate native solver processes.

At each EOS source break, close one-sided samples preserve finite changes of
equilibrium sound speed. A geometric enthalpy mesh is combined with the
uniform s mesh to resolve the low-sound-speed envelope.
All testing in this stage is *smooth EOS*: Maxwell sharp interfaces are
explicitly rejected, and no slow/rapid conversion condition is invented.

## Actual results: original TOV + Cowling versus native full path

All runs on the same connected Ubuntu 24.04 device with gfortran 13.3.0,
OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=1. `surface_h=1e-4`.
TOV: `rtol=1e-10, atol=1e-12, max_step_h=0.002`.
Cowling: `rtol=1e-9, atol=1e-12, ell=2`; matching roots at same
pre-identified brackets. Native source table `native_samples=8193`
before adding break / geometric knots.

| Source / central density (fm^-3) | M_native (Msun) | R_native (km) | Δf (Hz) | Δp1 (Hz) | Full path speedup |
|---|---:|---:|---:|---:|---:|
| primary, 0.35 | 0.7686502631 | 12.02829023 | -0.001708 | +0.080335 | 141.8x |
| primary, 0.53773398 | 1.4000000236 | 11.75207831 | -0.005033 | -0.008278 | 111.5x |
| primary, 0.70 | 1.7517946440 | 11.37843024 | -0.002856 | -0.026131 | 126.0x |
| extension medium, 1.15 | 2.0484564089 | 10.14174180 | -0.000340 | -0.000913 | 117.7x |
| extension fine, 1.15 | 2.0484563835 | 10.14174181 | -0.000920 | -0.009078 | 132.8x |

These ratios are the sum of one exact-Python TOV and two frequency-root
shootings versus one native-TOV build (including Python EOS source sampling),
native coefficient creation and two native root shootings, **after native
shared-library compilation was cached**. Brackets are pre-identified and
held fixed on both sides. The Python per-mode times come from the initial
exact-reference benchmark and are reused **only** because the Python reference
source and EOS hashes match; the final native numerics were independently
re-executed at the final code hash. This is NOT the first-install/cold-build
latency, and is not a full-native EOS generator.

All five cases have the same f and p1 W/V node counts as the exact reference.
Maximum |Δf|=0.005033 Hz, maximum |Δp1|=0.080335 Hz (the low-mass case).
Sampling convergence across 4097,8193,16385 shows persistent sensitivity of
the piecewise EOS p1 frequency. These residuals must not be interpreted
as exact-zero numerical errors. Per-case data, numerical parity residuals,
source hashes, run durations and convergence are in
`data/cowling_native_tov_20261119/FINAL_RECHECK.json`.

### Independent limits

- Constant-density acausal Schwarzschild interior (reference geometry only):
  the native TOV radius/mass agree with the exact solution to about 1e-13
  relative in the 4097–16385 sample tests.
- Analytic Gamma=2 polytrope (K=400, m_b=939):
  at n_c=0.4, independent Python DOP853 versus native f frequency differs by
  ~2.1e-6 Hz, and p1 by ~3.3e-6 Hz; both node counts agree. No inference
  about realism of this polytrope follows.
- TOV surface metric normalization, interior derivative constraints, and
  interchange of two native background objects are covered by tests.

## Scientific status and remaining work

**COWLING-NATIVE-02 status: SCOPED_IMPLEMENTATION_VERIFIED** for
Fortran TOV+Cowling integrators on the stated smooth EOS conditions,
the analytical limits and quoted numerical tolerances. Only
`numerically_checked`, not a blanket production qualification.
Full EOS source evaluation remains Python and its PCHIP representation
is a controlled approximation to original continuous functions.
On the low-mass piecewise EOS, default p1 precision is ~0.08 Hz; use
increased `native_samples` and report convergence for any tighter claims.

**Next selected research work:** `COWLING-MAXWELL-NATIVE-03`, native
sharp-Maxwell two-branch background/mode solver with separate **slow**
and **rapid** conversion junction conditions and independent interface
residual checks. Reference: Tonetto & Lugones,
Phys. Rev. D 101, 123029 (2020), DOI:10.1103/PhysRevD.101.123029.
This request does not authorize silently applying smooth formulas
through a density discontinuity.

## Reproduction

```bash
cd /home/cosmosapjw/ns_cowling_wt/cowling_native_tov_20261119
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=python_frontend
pytest -q python_frontend/tests/test_native_cowling_tov_20261119.py \
 python_frontend/tests/test_native_polytrope_cowling_20261119.py \
 python_frontend/tests/test_native_cowling_dvode_20261118.py \
 python_frontend/tests/test_stellar_background.py \
 python_frontend/tests/test_stellar_cowling.py \
 python_frontend/tests/test_cowling_reaudit_contracts.py \
 python_frontend/tests/test_css_cowling_adapter.py
# To recompute all expensive Python baseline solutions:
python3 scripts/benchmark_native_tov_cowling_20261119.py
# Recompute native at final source SHA against hash-locked Python baseline:
python3 scripts/recheck_native_tov_cowling_20261119.py
```

Unchanged tests/default APIs are not redefined as new physics evidence.
