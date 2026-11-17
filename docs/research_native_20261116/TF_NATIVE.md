# Native spatial Thomas–Fermi implementation

The Fortran modules `compose_eos_tf_models` and `compose_eos_spatial_tf` evaluate all eight original GPPVA bulk models, their finite-range scalar/vector fields, screened electrons, local density-dependent rearrangement, energy and chemical gradients. The constrained density/radius minimization and finite-family selection also run in Fortran. Python is retained as a test oracle and for launching/reporting the campaign.

`solve_tf_family_native(model_id,nbar,npoints,dimension,seed,selected,status,attempts)` returns every bounded minimization attempt. The first radii are 30/45/70 fm (45/70/100 below 0.003 fm^-3), with same-profile refinement when necessary. If strong clustered guesses all fail, a separately solved uniform native state supplies the mean species densities for explicit weak clustered guesses; a further minimization is mandatory. Such an initial state is not a substituted clustered EOS. The optional search is finite (two amplitudes, three radii and at most one refinement each). Failed initial/trial field solves retain `evaluation_status` and cannot resolve a family. The existing PDE, KKT, radius and search-bound criteria are unchanged.

At N16, all 8 models × 3 densities (0.003, 0.030, 0.075 fm^-3) × 6 geometry/seed families returned a stationary finite candidate: 144/144. The first campaign's 143/144 result is preserved separately. Its FSU2H 0.075 spherical-droplet initial guesses failed the vector Newton line search. The additional weak clustered initialization subsequently reached a clustered state with energy below its uniform initial state. This is an initialization failure and recovery, not evidence of a physical branch being absent.

The dedicated Python-launched native tests passed 27 cases, including independently differenced species and radius energy derivatives, all-eight homogeneous reductions, screened-electron PDE/neutrality/gauge checks, actual minimization for all models, and reevaluation of an archived FSU2R N384 profile. A separate compiled selector test passed four cases covering a finite minimum, a lower unconverged candidate, a missing family and a NaN energy.

These results establish native implementation and finite numerical samples. They do not certify a global minimum, continuum geometry transition, full density-support convergence, complete original author implementation, or an exact match to a supplied CompOSE table. The native API returns HOLD for unresolved families or feasible lower unconverged attempts.

Reproduction:

```sh
bash fortran_eos/tests/build_native_tf_test.sh /tmp/ns-native-tf-test
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /tmp/ns-native-tf-test/native_tf_selection_test
PYTHONPATH=python_frontend OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q python_frontend/tests/test_native_spatial_tf.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 data/research_native_20261116/tf/run_family_campaign.py --driver /tmp/ns-native-tf-test/native_tf_driver --output /tmp/ns-native-tf-new-campaign
```

`data/research_native_20261116/tf/RECEIPT.json` records source and evidence hashes. Per-family stdout stores the selected profile, stderr stores all attempt summaries, and JSON stores numeric results and execution identity. The production caller receives the full native attempt states, including profiles where available.
