# NS_Cowling p1 backend selector: physical admission regression

Date: 2026-11-20 KST
Parent: `f96936bb41f200de3d43b357da4506f4ddf93424`
Branch: `research/cowling-mode-selector-validation-20261120`
Scope: C6/C5 bounded validation and repair. Numerical ODE and EOS physics unchanged.

## Reproduced defect
`scripts/summarize_hybrid_accuracy_20261120.py::choose_verified_backend`
formerly trusted a copied `MATCHED_REFERENCE_PASS` field and accepted
the fastest candidate whenever its stored frequency difference and relative
join were below thresholds. It did not inspect the required p1 mode nodes,
negative join data, or a frequency-difference inconsistent with its source
reference. A fabricated fast invalid candidate won against a valid one.

Before repair, the new regression suite observed four expected assertion
failures and one passing control (`4 failed, 1 passed`).

## Implemented correction
The selector now requires:
* explicit p1 `nodes == [1,1]` for both recorded p1 benchmark cases;
* finite positive measured time and frequencies;
* finite relative junction residual in [0, 1e-5];
* a recorded signed frequency difference equal to (f_native - f_reference)
  within subtraction roundoff (max(1e-9 Hz, 32 binary ulps));
* finite signed error within the caller's positive error threshold.

An invalid candidate is skipped, and the original no-candidate state is
retained when all candidates fail. No global mode-identity guarantee or
independent error bound is inferred from the nodes alone.

## Fresh checks
From this worktree, with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`
and `PYTHONPATH=python_frontend`:

```sh
python3 -m pytest -q \
 python_frontend/tests/test_cowling_empirical_selector_physics_20261120.py \
 python_frontend/tests/test_sharp_cowling_hybrid_cs2_20261120.py \
 python_frontend/tests/test_sharp_cowling_hybrid_evidence_20261120.py \
 --run-slow --disable-warnings
```

Observed: `14 passed in 25.08s`; both physical M3 slow / M2 rapid
hybrid tests included. The prior fixed benchmark JSON was read without
mutation; all six selection decisions at 0.05, 0.005 and 0.001 Hz
exactly matched the previously recorded decisions. Python compileall
and `git diff --check` succeeded.

## Independent interface-law check (no code modification)
With `w=W/r^(ell+1)`, `v=V/r^ell`,
`a=omega_geom^2 exp(lambda-2Phi) r`, `b=dPhi/dr`,
and `H=epsilon+P_t`, the common traction has the form
`H(a*v + b*w)`. Imposing its continuity together with `[w]=0`
for slow conversion or `[v]=0` for rapid conversion yields exactly
the transfer maps implemented in
`fortran_eos/src/compose_eos_cowling_dvode.f90`.
The symbolic equalities and zero interface residuals were evaluated
using Wolfram Language. This is algebraic verification of the encoded
interface mapping, not a physical validation of other approximations.
Primary reference: Tonetto and Lugones, Phys. Rev. D 101, 123029 (2020),
https://arxiv.org/abs/2003.01259

## Identity and scope
The old `data/cowling_hybrid_geometry_20261120/verification/VERIFICATION_RECEIPT.json`
attests to the immutable parent revision, not this modified selector.
Historical `RESULT.json`, `DECISIONS.json` and the source-locked EOS
bundle were not rewritten. This bounded new regression does not claim
a complete 203-test rerun, broader EOS/mass/conversion coverage, or a
global frequency error bound.
