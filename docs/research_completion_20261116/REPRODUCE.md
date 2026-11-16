# Reproduce the independent EOS components

Run from the repository root in a source checkout. Python3.12, NumPy2.4.2,
SciPy1.17.0 and gfortran13.3.0 were used locally. Install the test dependencies
with `python3 -m pip install -e '.[test]'`; mpmath supplies independent
high-precision tests. These commands calculate physical states; they do not
copy the target CompOSE thermodynamics.

```sh
bash fortran_eos/build.sh build/compose_eos_fortran
export PYTHONPATH=python_frontend
export COMPOSE_EOS_FORTRAN_ROOT=.
export COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
```

Use fresh output directories so old results cannot masquerade as a new run.

| Component | Command | Meaning |
|---|---|---|
|All five cores|`python3 scripts/generate_core_eos.py build/replay_core --points 161`|Caller-defined grids; composition and native readback. BBB uses its declared published-energy reconstruction.|
|All five SPG candidates|`python3 scripts/generate_spg_crust.py build/replay_spg --points 513`|Spherical CLD, explicit candidate surface calibration, complete gas proton/neutron variables and grid comparison.|
|Outer crust|`python3 scripts/build_outer_crust.py build/replay_outer --points 513 --coarse-points 257`|Independent mass-catalogue Gibbs minimization and transition/drip roots.|
|BSk22 lower inner crust|`python3 scripts/build_bsk22_lower_inner.py build/replay_bsk22 --points 129`|Published analytic EOS/composition; separate source-pressure and energy-derivative modes.|
|Optional BSk22 boundary reconstruction|`python3 scripts/build_bsk22_reconstruction.py build/replay_bsk22_anchored`|Preserves source pressure and integrates the first law from the independently computed outer boundary. Approximate composition remains labeled.|
|All eight spatial TF models|`python3 scripts/build_gppva_tf.py build/replay_tf --workers 6`|Finite-volume meson/Coulomb fields, screened electrons, geometry and radius optimization.|
|TF mesh refinement|`python3 scripts/refine_gppva_tf.py build/replay_tf build/replay_tf_refined --workers 6`|Actual profile reoptimization on finer meshes; retain unresolved states explicitly.|

The actual dense TF campaign additionally supplied:

```sh
python3 scripts/build_gppva_tf.py build/replay_tf_dense --workers 6 \
  --resolutions 48 96 \
  --densities .0002 .00035 .0005 .00075 .001 .0015 .002 .0035 .005 .0075 \
  .01 .015 .02 .025 .03 .035 .04 .045 .05 .055 .06 .065 .07 .075 .08 .085 .09 .095 .1
python3 scripts/refine_gppva_tf.py build/replay_tf_dense build/replay_tf_dense_refined --workers 6
```

An optional `--warm-start PREVIOUS_CAMPAIGN` changes only TF initial guesses;
the original trajectory and all source hashes remain recorded in each receipt.
A fresh run may choose a different local stationary branch. Compare energy,
constraints and numerical errors rather than requiring optimizer byte identity.

The TF native table contains electron fractions and honest availability masks.
Read `tf_composition.py` companion profiles for spatial neutron/proton/electron
densities and mean species fractions. Do not interpret absent discrete nuclear
A/Z or bound/free decomposition as zero nuclei. Do not interpolate across a
density that the campaign explicitly held.

The BSk22 author-code comparison needs the external original Fortran file;
its URL, SHA256, driver and compilation instructions are retained in
`data/research_completion_20261116/bsk22_lower_inner/SOURCE.json`.
The source has no explicit redistribution license in its header, so it is not
vendored. The Python equations and coefficients are independently implemented.
The source mass inventory required by outer crust is committed under
`data/outer_crust/`, with acquisition identities and the replay intake script.
The checkout layout is currently required; a non-editable wheel does not yet
package every research input file.

## Verification

```sh
PYTHONPATH=python_frontend python3 -m pytest -q python_frontend/tests pipeline/tests
PYTHONPATH=eos_generator:eos_generator/compose_eos_pkg python3 -m pytest -q eos_generator/tests
python3 data/research_completion_20261116/core/check_causal_boundaries.py
sha256sum -c SHA256SUMS_v4.txt
```

The core entry point is also executed in CI and its newly generated tables are
uploaded. The expensive multi-start TF campaigns are retained with replayable
scripts and independent reviews; ordinary unit tests do not replace them.

Outer/inner matching prescriptions and candidate rest-mass reconstruction are
explicit research variants, documented in the physics review. The commands above
do not silently turn independently generated components into a combined
core–crust EOS. No new stellar calculation is part of this loop.
