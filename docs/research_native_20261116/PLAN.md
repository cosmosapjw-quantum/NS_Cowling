# Fortran EOS generation execution plan

Baseline: 4de49430560930be84d9a4147d71c3c9841da2db. The current user request authorizes theory and coding toward original core5/crust13 regeneration, with Fortran physics and Python frontend only. Existing dirty checkouts remain untouched.

Architecture: native model functionals, gradients, fields, minimization, phase selection and physical acceptance; Python transports input, launches binaries, validates file identities and presents results. Existing Python solvers remain archival/test oracles. No EOS table replay or proxy-model substitution.

- [x] Native core5: APR two-charge Gibbs, BBB published-energy reconstruction, GM1, KDE0v1, SLy4.
- [x] Native GDFM/SPG5: explicit surface inputs, gas-proton KKT, drip/endpoint and EOS states.
- [x] Native spatial TF8: model-specific couplings, finite-volume fields, neutrality, geometry/seed search, morphology and profile diagnostics.
- [x] Pinned licensed native SLSQP service and constrained optimizer tests.
- [x] Native nuclear-mass outer crust and explicitly approximate BSk analytic lower crust.
- [x] Primary-source assessment: lower attachment resolved as a documented boundary; original SPG fit explicitly unresolved.
- [x] Native build/CLI, thin frontend, actual full-cohort execution, source/checksum receipts, bounded independent verification, publication.

Each implementation owns disjoint modules/tests and records failing-before/working-after tests. Source fidelity, physical identities, numerical support and CompOSE comparison are independent statuses. Missing original fit stays unresolved; explicit candidates are executable but not original-equivalent. No full BHF or ETFSI claim from energy tables or analytic fits. No arbitrary relative energy offsets, hidden failed densities, global-minimum claims from finite searches, or new stellar/core-crust campaign.

Reuse prior156-node and subsequent NL3 evidence when identity permits. New native code requires actual native numerical tests and regeneration; changing the producer invalidates Python execution claims as evidence of native execution. Parent integrates shared build/frontend and pushes only after verification; no force push or merge.
