#!/usr/bin/env bash
# Compile only the actual TF diagnostic producers, in an isolated directory.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD="${1:-$ROOT/build/blocker_tf_probes}"
mkdir -p "$BUILD"
FC="${FC:-gfortran}"
"$FC" -std=legacy -O2 -ffixed-line-length-none -c "$ROOT/fortran_eos/vendor/slsqp/slsqp_optmz.f" -o "$BUILD/slsqp.o"
for module in compose_eos_types compose_eos_optimize_native compose_eos_tf_models compose_eos_spatial_tf; do
 "$FC" -std=f2018 -O2 -ffree-line-length-none -fcheck=all -J"$BUILD" -I"$BUILD" -c "$ROOT/fortran_eos/src/$module.f90" -o "$BUILD/$module.o"
done
COMMON=("$BUILD/compose_eos_types.o" "$BUILD/compose_eos_optimize_native.o" "$BUILD/compose_eos_tf_models.o" "$BUILD/compose_eos_spatial_tf.o" "$BUILD/slsqp.o")
for name in blocker_tf_branch_derivative blocker_tf_common_tangent blocker_tf_common_tangent_polish; do
 "$FC" -std=f2018 -O2 -ffree-line-length-none -fcheck=all -I"$BUILD" "$ROOT/fortran_eos/tests/$name.f90" "${COMMON[@]}" -llapack -lblas -o "$BUILD/$name"
done
"$FC" -std=f2018 -O2 -ffree-line-length-none -fcheck=all -I"$BUILD" "$ROOT/fortran_eos/tests/blocker_tf_export_chord.f90" "$BUILD/compose_eos_types.o" -o "$BUILD/blocker_tf_export_chord"
