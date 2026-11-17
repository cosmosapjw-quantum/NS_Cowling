#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD="${1:-$ROOT/build/native_tf_test}"
FC="${FC:-gfortran}"
mkdir -p "$BUILD"
"$FC" -std=legacy -O2 -ffixed-line-length-none -c "$ROOT/fortran_eos/vendor/slsqp/slsqp_optmz.f" -o "$BUILD/slsqp.o"
COMMON=("$ROOT/fortran_eos/src/compose_eos_types.f90" "$ROOT/fortran_eos/src/compose_eos_optimize_native.f90" "$ROOT/fortran_eos/src/compose_eos_tf_models.f90" "$ROOT/fortran_eos/src/compose_eos_spatial_tf.f90")
for PROGRAM in native_tf_driver native_tf_selection_test; do
 "$FC" -std=f2018 -O2 -ffree-line-length-none -fcheck=all -J"$BUILD" -I"$BUILD" "${COMMON[@]}" "$ROOT/fortran_eos/tests/$PROGRAM.f90" "$BUILD/slsqp.o" -llapack -lblas -o "$BUILD/$PROGRAM"
done
printf '%s\n' "$BUILD/native_tf_driver" "$BUILD/native_tf_selection_test"
