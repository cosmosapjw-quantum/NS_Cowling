#!/usr/bin/env bash
# Build the finite coexistence diagnostics without changing production kernels.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD="${1:?build directory required}"
mkdir -p "$BUILD"
BUILD="$(cd "$BUILD" && pwd)"
bash "$ROOT/fortran_eos/tests/blocker_tf_build.sh" "$BUILD"
for name in tf_finish_diagnostic tf_finish_chord tf_finish_single_state tf_finish_profile_readback; do
 gfortran -std=f2018 -O2 -ffree-line-length-none -fcheck=all -I"$BUILD" \
  "$ROOT/fortran_eos/tests/$name.f90" "$BUILD/compose_eos_types.o" \
  "$BUILD/compose_eos_optimize_native.o" "$BUILD/compose_eos_tf_models.o" \
  "$BUILD/compose_eos_spatial_tf.o" "$BUILD/slsqp.o" -llapack -lblas -o "$BUILD/$name"
done
printf '%s\n' "$BUILD/tf_finish_diagnostic" "$BUILD/tf_finish_chord" "$BUILD/tf_finish_profile_readback"
