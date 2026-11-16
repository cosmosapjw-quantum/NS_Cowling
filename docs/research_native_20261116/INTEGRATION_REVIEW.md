# Native generation integration review

Base: `4de49430560930be84d9a4147d71c3c9841da2db`. Review date: 2026-11-16. Scope: the new Fortran generator entry point, public Python dispatcher, native readback/provenance transport and boundary representation. This is a bounded implementation review, not a new proof of the nuclear functionals or global cell minima.

## Findings and repairs

1. **Uniform TF export:** the first app version marked the complete charge inventory while leaving Abar=0. The native writer correctly refuses that combination. The app now uses the existing empty-cluster carrier Abar=1, Zbar=0, Yh=0 for uniform matter. Nonuniform TF free/cluster inventories remain masked as unavailable rather than being inferred from total nucleon fractions.
2. **TF mass headers:** the first app version used 939 MeV universally. It now obtains both mass headers from the actual native model, including unequal DD2 nucleon masses and the 938 MeV TM1e convention. The thermodynamic quantities and the header normalization therefore share the declared model convention.
3. **Failure provenance:** attempt CSV now includes the escaped native evaluation status, optimizer mode, iteration count and evaluation count. Failed profile trials remain recorded; unresolved families or a lower feasible unconverged trial prevent a completed-table status.
4. **Writer checks:** blanket `allow_source_physical_violations=.true.` was removed from the generator. Fresh native rows pass the writer's ordinary composition-preservation and applicable closure checks.
5. **Receipt trust:** an independent reproducer showed that nonempty files containing only “not a native table” could be accepted by the first Python verifier if accompanied by a success-shaped RESULT.json. The revised verifier invokes the actual Fortran reader and records native C1–C11 diagnostics. It checks the model/family/backend identity, actual row count, requested core grid (allowing APR's inserted in-window endpoints), every requested TF grid and mesh directory, and equality of the main TF table to the finest requested mesh. Lower segments receive their own native readback and hashes. This repair validates transport and completeness; Python does not evaluate a functional or select a phase.

The deliberate red reproducer is retained in `data/research_native_20261116/integration/independent_verifier_red.json`. Its acceptance was a verifier defect, not a successful physical calculation.

## Observed verification

- Twelve dedicated transport tests pass against the real Fortran reader and diagnostics, including malformed tables with a fabricated receipt, wrong density with unchanged row count, an omitted requested mesh, disagreement with the finest mesh and malformed lower-crust composition.
- The real public command `python -m compose_eos.cli generate GM_GM1 … --points 3 --n-min .1 --n-max .4` completes with `REQUESTED_ROWS_NATIVE_VERIFIED`, return code 0 and three native rows. Both generator readback and independent frontend-triggered Fortran readback succeed.
- Logs and the GM1 receipt are under `data/research_native_20261116/integration/`. A prior attempted `python -m compose_eos` invocation failed because the package has no `__main__`; that command error is preserved separately and was corrected using the existing CLI module.
- The legacy linker warning about an executable stack from native callback closures is preserved in the execution log; this review did not change that implementation.

## Physical interpretation and remaining limits

Native consistency results are diagnostics, not a requirement to numerically reproduce a CompOSE target or to make every finite-difference test pass on sparse support. No Python physical-domain calculation was introduced into the production path. A successfully regenerated source representation can still have limited physical support, acausality, sparse derivative support or a discontinuous seam; those findings must be interpreted separately.

Outer, BSk22 lower and model-specific inner data are exported as distinct segments. The tiny outer/BSk density overlap and pressure/chemical-potential differences are retained, with cross-segment interpolation explicitly unapproved. The documented nominal 0.002 fm^-3 attachment is not silently replaced by a newly imposed common tangent. Finite TF samples and phase changes do not establish a continuum mesh limit, a continuous transition density or a global minimum.

SPG surface coefficients remain explicit candidate inputs; regenerability does not establish equivalence to the original author calibration. The parent campaign supplies the all-model numerical evidence and final publication status. This review alone does not claim that all original core/crust scientific equivalence questions are closed.
