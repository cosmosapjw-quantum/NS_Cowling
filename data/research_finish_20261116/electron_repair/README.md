# Electron reduction evidence

Final result: `RESULT.json`; derivation: `docs/research_finish_20261116/ELECTRON_REPAIR.md` from repository root. Final production SHA-256 is `8be3c6756f1b04a38cb59fdeeacb5c0aed5c98d8337c8ed14c64e7091501bdca`.

- `baseline_spatial_tf_equilibrium.py`, `RED_RECEIPT.json`, `red.log`, `red_tests.log`, `warm_failure_input.json`, `cold_failure_input.json`: original failed DD2 implementation, two full minimizations, exact failed electron inputs and two red regression tests. The saved source is imported by `reproduce_failure.py`, independently of current production.
- `neutral_pinned_spatial_tf_equilibrium.py`, `SLAB_RED_RECEIPT.json`, `slab_red.log`, `slab_red_tests.log`, `GPPVA_*_slab_failure_input.json`: first repair and newly exposed six slab failures. These six records contain five distinct electron inputs; the cold inputs coincide. `reproduce_slab_failures.py` imports the saved first repair.
- `dilute_stopping_red.log`: two counterexamples motivating a tighter near-vacuum internal target.
- `GREEN_RECEIPT.json`, `REPAIRED_CELLS.json.gz`, `green_cells.log`: final-source DD2 rod replays with full profiles and all attempts.
- `FINAL_SLAB_CELLS.json.gz`, `final_slab_cells.log`: final-source six TW/DDME2 slab minimizations, all stationary under the existing criteria.
- `green_tests.log`: final 53-test result, including full-functional derivative checks and all eight actual failure records.
- `INITIAL_*`, `PINNED_*`, `initial_*`, `green_initial.log`: retained intermediate evidence, not final-source claims.
- `FLUX_*`, `flux_*`, `SLAB_FLUX_CELLS.json.gz`, `slab_flux_cells.log`, `prototype_flux_repair.py`, `replay_flux_cells.py`, `check_flux_meshes.py`: conservative-flux prototype, successful actual slab replays, and a recorded N768 rounding limitation.
- `EXTENDED_FLUX_MESH_CHECKS.json`, `extended_*`, `check_extended_flux_meshes.py`: bounded extended-potential-precision prototype resolving all18 N192/N384/N768 interpolated-profile checks. Final production adds explanatory comments to this tested algorithm.

Replays from repository root (these scripts write new run outputs):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 data/research_finish_20261116/electron_repair/reproduce_failure.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 data/research_finish_20261116/electron_repair/reproduce_slab_failures.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 data/research_finish_20261116/electron_repair/replay_repaired_cells.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 data/research_finish_20261116/electron_repair/replay_final_slab_cells.py
```

Failed input arrays and seed archive hashes are retained; final profiles are finite JSON data. Source identities are exact byte identities. Floating optimization results can vary across environments and gzip container bytes include write metadata; physical residual criteria and seed/source identities are the scientific reproducibility contract, not identical compressed-output bytes. This repair neither calibrates a nuclear model nor establishes an EOS-wide or continuum result.
