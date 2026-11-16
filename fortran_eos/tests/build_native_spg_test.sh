#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/../.." && pwd)"
out="$1"
mkdir -p "$out"
gfortran -std=legacy -O2 -ffree-line-length-none -ffixed-line-length-none -fcheck=all -J"$out" -I"$out" \
 "$root/fortran_eos/src/compose_eos_types.f90" \
 "$root/fortran_eos/vendor/slsqp/slsqp_optmz.f" \
 "$root/fortran_eos/src/compose_eos_optimize_native.f90" \
 "$root/fortran_eos/src/compose_eos_gdfm_native.f90" \
 "$root/fortran_eos/src/compose_eos_spg_native.f90" \
 "$root/fortran_eos/tests/native_spg_probe.f90" -o "$out/native_spg_probe"
