# Executed spatial TF campaign

All eight models have 23/23 selected stationary source-interval nodes from 0.002 to 0.1 fm^-3. The final tables contain 226 rows including lower-density diagnostic extensions. Uniform phases at the upper end are explicitly identified. No held density is interpolated into the table.

| Model | Final rows | Lower-extension held densities (fm^-3) | Last clustered / first uniform sampled n | Max source-interval mesh P change |
|---|---:|---|---|---:|
| GPPVA_DD2 | 28 | 0.0002 | 0.07 / 0.075 | 0.1922% |
| GPPVA_DDME2 | 27 | 0.0002, 0.00035 | 0.07 / 0.075 | 0.1729% |
| GPPVA_TW | 28 | 0.0015 | 0.075 / 0.08 | 0.1724% |
| GPPVA_FSU2 | 28 | 0.0002 | 0.05 / 0.055 | 0.1780% |
| GPPVA_FSU2R | 29 | none | 0.09 / 0.095 | 0.3385% |
| GPPVA_FSU2H | 29 | none | 0.095 / 0.1 | 0.1704% |
| GPPVA_NL3wrL55 | 29 | none | 0.08 / 0.085 | 0.1000% |
| GPPVA_TM1e | 28 | 0.00075 | 0.095 / 0.1 | 0.2101% |

The transition entries are sampled geometry brackets, not a high-precision phase-transition solve. Six starting geometry classes were compared at N48/N96; N192 refined every observed geometry within 0.01 MeV/baryon of the coarse minimum plus uniform matter. N384 was used below and at 0.002 fm^-3. This is a finite search, not a proof of the global minimum.

The radial-mesh statement is independent of density interpolation accuracy. Final SUPPORT.json files retain actual midpoint pressure/secant diagnostics; the independent tf_density_review collection assesses a fixed-N192 density grid with pressure, energy and electron-fraction errors and phase labels.

Retained stages:

- campaign_v1: 1,440 initial three-species optimizations; weakly weighted electron residual failures are explicit.
- campaign_v2: 2,784 attempted electron-eliminated cells; stale electron warm-start failures are explicit.
- campaign_v3: 2,784 attempted cells after neutral warm-start repair, before final bracketing padding.
- campaign_refined_initial: unchanged N192/N384 profiles before numerical stopping/bound-seed repairs.
- final: native tables, total composition companions, 226 selected profiles and all source-interval repair trials.
- source_versions: SHA-verified prior numerical implementations and original refinement runner. Current physics modules are tracked in the repository.
- logs: original failures, repaired runs and test logs. replay_adversarial.py reproduces the stale-electron and endpoint-bracket failures, then checks the repaired solver.
- numerical_counterexamples.json.gz: actual failed energy-floor and nearly absent-species profiles used by two regression tests.

The 32 original/expanded TF tests passed. Two additional replay regression tests passed after introducing the refinement restart script. The parent full-suite run is recorded separately; these targeted counts are not relabeled as a combined-suite run.

Independent reevaluation of all final 226 profiles is in ../gppva_bulk/screened_campaign_review.json. Original 222-profile review identities remain in screened_refined_campaign_review.json. The all8 functional and optimized density-pressure derivations are in ../tf_review/.

Reproduce the main route from the repository root with a fresh output directory:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python3 scripts/build_gppva_tf.py OUTPUT --workers 6 --densities .0002 .00035 .0005 .00075 .001 .0015 .002 .0035 .005 .0075 .01 .015 .02 .025 .03 .035 .04 .045 .05 .055 .06 .065 .07 .075 .08 .085 .09 .095 .10
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python3 scripts/refine_gppva_tf.py OUTPUT REFINED_OUTPUT --workers 6
```

The current refinement runner already includes the physical-criteria-preserving restarts learned from this campaign. finalize_gppva_tf.py records how the historical initial refinement was repaired and packaged; it is not an independent physical model.
