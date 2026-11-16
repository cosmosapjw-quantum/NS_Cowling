# Isolated coupled-vector band solve

`PROPOSED.patch` is the production proposal; it was **not applied by this experiment**. `baseline_spatial_tf.py` and `band_spatial_tf.py` preserve the exact compared sources. `coupled_band_function.txt` is the helper included in the patch. The corresponding electron modules are also preserved so the full isolated calculation can be replayed.

`BENCHMARK.json` and `benchmark.log` report 12 paired evaluations on retained FSU2R/FSU2H/TM1e N192 profiles and explicitly labelled N384 resamples, with cold/perturbed field starts and five timing samples. They also retain five singular/pivot/indefinite/matrix-structure adversaries. `FSU2R_N384_CELL.json.gz` is an actual new N384 minimization; `FSU2R_N384_RECEIPT.json` retains its source/input hashes, physical criteria, timing and dense/band final-profile reevaluation. `prototype_tests.log` is the 63-test isolated-prototype result. `missing_helper_red.log` only demonstrates that the proposed helper does not exist in baseline production; it is not a failed baseline physics check.

Reproduce from repository root (scripts write new outputs):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 data/research_finish_20261116/vector_band_prototype/benchmark_vector_fields.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 data/research_finish_20261116/vector_band_prototype/run_fsu2r_cell.py
```

Derivation and scope: `docs/research_finish_20261116/VECTOR_BAND_PROTOTYPE.md`. New focused tests: `python_frontend/tests/test_spatial_tf_vector_band.py`. The linear operation is algebraically equivalent; no physics, stationarity criterion or source-table gate changes. Timing ranges describe this loaded machine and do not predict whole-campaign speedup. SHA256SUMS identifies evidence bytes; gzip output bytes may change on replay without changing scientific content.
