#!/usr/bin/env bash
# Standalone Fortran C12 scientific test build from exact checked-out sources.
# Modified project modules use -fcheck=all. Frozen Netlib F77 receives
# compatibility flags, without a false claim of vendor-wide memory safety.
set -euo pipefail
if [[ $# -ne 1 ]];then
  echo "usage: $0 /abs/path/to/c12_fortran_driver" >&2;exit 2
fi
OUT=$1
[[ "$OUT" = /* ]] || exit 2
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
DIR=$(dirname "$OUT")
mkdir -p "$DIR"
cd "$DIR"
gfortran -O0 -g -fcheck=all -ffree-line-length-none -c \
 "$ROOT/fortran_eos/src/compose_eos_cowling_dvode.f90" \
 "$ROOT/fortran_eos/src/compose_eos_cowling_source_exact.f90" \
 "$ROOT/fortran_eos/src/compose_eos_cowling_sharp_source_root.f90" \
 "$ROOT/fortran_eos/tests/cowling_source_exact_sharp_root_c12.f90"
gfortran -O2 -std=legacy -c \
 "$ROOT/fortran_eos/vendor/odepack/vode.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgefa.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgesl.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgbfa.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgbsl.f"
gfortran -g compose_eos_cowling_dvode.o \
 compose_eos_cowling_source_exact.o \
 compose_eos_cowling_sharp_source_root.o \
 cowling_source_exact_sharp_root_c12.o \
 vode.o dgefa.o dgesl.o dgbfa.o dgbsl.o -lblas -o "$OUT"
sha256sum "$OUT"
