#!/usr/bin/env bash
# Supply a new output directory; retained source scripts stay in the checkout.
set -euo pipefail
export NS_TF_CONTINUATION_OUTPUT="${1:?new output directory required}"
if [ -e "$NS_TF_CONTINUATION_OUTPUT" ]; then echo "Output already exists" >&2; exit 2; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="python_frontend:/tmp/ns_cowling_audit_20261117_pydeps"
PYTHON="/home/cosmosapjw/ns_cowling_wt/.venv-gen/bin/python"
B=build/tf_continuation_20261117
D=data/tf_pair_expansion_20261117/tf_continuation
bash fortran_eos/tests/blocker_tf_build.sh "$B"
gfortran -std=f2018 -O2 -ffree-line-length-none -fcheck=all -I"$B" fortran_eos/tests/tf_continuation_20261117.f90 "$B/compose_eos_types.o" "$B/compose_eos_optimize_native.o" "$B/compose_eos_tf_models.o" "$B/compose_eos_spatial_tf.o" "$B/slsqp.o" -llapack -lblas -o "$B/tf_continuation_20261117"
"$PYTHON" scripts/continue_tf_20261117.py --phase coarse
"$PYTHON" scripts/continue_tf_20261117.py --phase refine
"$PYTHON" "$D/repair_run.py"
"$PYTHON" "$D/fine_run.py"
"$PYTHON" "$D/endpoint_run.py"
"$PYTHON" "$D/derivative_run.py"
"$PYTHON" scripts/continue_tf_20261117.py --phase validate
"$PYTHON" "$D/summarize.py"
"$PYTHON" -m pytest -q python_frontend/tests/test_tf_continuation_20261117.py
