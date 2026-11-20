#!/usr/bin/env bash
# Rebuild standalone sharp one-Maxwell native regression on this source tree.
# Modified project Fortran is fully bounds-checked; frozen Netlib F77
# compatibility units are compiled separately without incompatible bounds
# instrumentation. This does NOT certify vendor source memory safety.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "usage: $0 /absolute/fortran-test-binary" >&2;exit 2
fi
OUT=$1
[[ "$OUT" = /* ]] || exit 2
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
DIR=$(dirname "$OUT")
mkdir -p "$DIR"
cd "$DIR"
gfortran -O0 -g -fcheck=all -ffree-line-length-none -c \
  "$ROOT/fortran_eos/src/compose_eos_cowling_dvode.f90" \
  "$ROOT/fortran_eos/src/compose_eos_cowling_sharp_root.f90" \
  "$ROOT/fortran_eos/tests/cowling_sharp_onejump_root_test.f90"
gfortran -O2 -std=legacy -c \
  "$ROOT/fortran_eos/vendor/odepack/vode.f" \
  "$ROOT/fortran_eos/vendor/odepack/dgefa.f" \
  "$ROOT/fortran_eos/vendor/odepack/dgesl.f" \
  "$ROOT/fortran_eos/vendor/odepack/dgbfa.f" \
  "$ROOT/fortran_eos/vendor/odepack/dgbsl.f"
gfortran -g compose_eos_cowling_dvode.o \
  compose_eos_cowling_sharp_root.o cowling_sharp_onejump_root_test.o \
  vode.o dgefa.o dgesl.o dgbfa.o dgbsl.o -lblas -o "$OUT"
sha256sum "$OUT"
