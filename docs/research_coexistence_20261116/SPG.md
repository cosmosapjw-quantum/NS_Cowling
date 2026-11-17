# SPG source-row surface-derivative check

Baseline HEAD: `a0b158f3dc4af9825192bf6a185a6331619868a5`.

## Result

A new native Fortran finite-difference calculation evaluated an explicit
SPG CLD functional using the inferred surface parameters at the 760 eligible
source-composition states. At fixed `(n_b, f, n_i, R)`, it perturbed only
the dense-phase proton fraction
`y=Z/A`. This is the partial variation relevant to the beta stationarity
condition. It directly included native GDFM bulk energy, the native uniform
electron energy, cell Coulomb energy, and the inferred Scurto surface form.

For each state, centered steps `1e-3, 3e-4, 1e-4, 3e-5, 1e-5` were used.
Across all 760 rows, the RMS difference between the finite-difference beta
derivative and its analytic full-derivative expression decreased from
`2.5183e-4 MeV` at step `1e-3` to `5.6690e-8 MeV` at `1e-5`; the maximum
at the smallest step was `1.9543e-7 MeV`. The trend is consistent with
second-order centered-difference convergence followed by floating-point
roundoff. This independently checks the analytic derivative implementation
against the same explicit energy functional.

Using the diagnostic replacement of the surface contribution by
`(1-f) sigma'(y)` reduces the RMS table beta residual from
`0.32077 MeV` to `0.004496 MeV` and its maximum from `0.78351 MeV` to
`0.008099 MeV`. This quantifies the already observed pattern. It does not
establish that the original authors used this prescription, that their code
omitted a term, or that the remaining residual has a unique cause. The fit
coefficients remain inferred from source composition under the stated
full-`A`, proton-free-gas convention; their identity with hidden source fit
inputs is unproved.

No production physics, default coefficients, or source tables changed.
This is a scoped derivative check, not source-equilibrium equivalence or an
independent EOS validation.

## Reproduction and identity

From the repository root:

```sh
python3 fortran_eos/tests/spg_residual_run.py /tmp/spg_native_fd_reproduction
```

The wrapper first rechecks all 15 pinned source-table hashes, then writes the
760-row native input, build/run stdout and stderr, finite-difference CSV, and
`SUMMARY.json` in the specified fresh output directory. It compiles
`spg_residual_fd.f90` against
`compose_eos_types.f90` and `compose_eos_gdfm_native.f90` using
`gfortran -std=legacy -O2 -ffree-line-length-none -fcheck=all`, then runs
the native program with OMP/OpenBLAS/MKL thread counts set to 1.

The native driver checks every model load and each GDFM baryon evaluation
before using its output. A deliberate invalid-density probe returned status
21 and recorded `FAIL dense reference SPG_M1 999999 invalid GDFM input`,
confirming that this error path exits nonzero. The first successful raw CSV
from before that guard was added remains at
`data/research_coexistence_20261116/spg/native_fd_reproduction/stdout.csv`;
its SHA256 is the same as the final guarded run.

The exact source-table hashes for all 15 SPG `eos.nb`, `eos.compo`, and
`eos.thermo` files are copied into `SUMMARY.json` and were rechecked:
15/15 matched. The prior inverse-row input, output, and inferred-parameter
hashes, direct producer/helper source hashes, native executable hash, input
manifest hash, and result CSV hash are also in that receipt. The recorded
result CSV SHA256 is
`b9ccf639b61150f238cc8b90d8034aa0d6853b8798c2a0598b986d28921f8160`.

Two setup attempts stopped before native execution because they included
120 rounding-probe records in the expected 760-row inventory. Those attempts
are preserved in `FAILED_ATTEMPTS.log`; the corrected wrapper selects the
760 original source rows and joins them by model and row identifier. There
were no native execution failures in the completed run.
