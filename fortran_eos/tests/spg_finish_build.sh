#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/../.." && pwd)"
build="${1:?isolated native build directory required}"
mkdir -p "$build"
build="$(cd "$build" && pwd)"
bash "$root/fortran_eos/build.sh" --ensure "$build"
for name in spg_finish_validation spg_finish_regression spg_finish_branch_probe spg_finish_drip_curvature; do
 gfortran -std=f2018 -O2 -ffree-line-length-none -fcheck=all -I"$build/mod" -J"$build/mod" \
  "$root/fortran_eos/tests/$name.f90" \
  "$build/obj/compose_eos_types.o" "$build/obj/compose_eos_gdfm_native.o" \
  "$build/obj/compose_eos_optimize_native.o" "$build/obj/compose_eos_spg_native.o" \
  "$build/obj/slsqp_optmz.o" -o "$build/$name"
done
