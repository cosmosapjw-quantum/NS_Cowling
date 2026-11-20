#!/usr/bin/env bash
# Build independent direct-source EOS test from the current checked-out code.
# Project Fortran receives -fcheck=all. Legacy frozen Netlib F77
# wrappers use -std=legacy without incompatible whole-vendor bounds checks.
set -euo pipefail
if [[ $# -ne 1 ]]; then echo "usage: $0 /abs/path/to/test-executable" >&2;exit 2;fi
OUT=$1
[[ "$OUT" = /* ]] || exit 2
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
DIR=$(dirname "$OUT")
mkdir -p "$DIR"
cd "$DIR"
gfortran -O0 -g -fcheck=all -ffree-line-length-none -c \
 "$ROOT/fortran_eos/src/compose_eos_cowling_dvode.f90" \
 "$ROOT/fortran_eos/src/compose_eos_cowling_source_exact.f90" \
 "$ROOT/fortran_eos/tests/cowling_exact_source_cs2_test.f90"
gfortran -O2 -std=legacy -c \
 "$ROOT/fortran_eos/vendor/odepack/vode.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgefa.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgesl.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgbfa.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgbsl.f"
gfortran -g compose_eos_cowling_dvode.o \
 compose_eos_cowling_source_exact.o cowling_exact_source_cs2_test.o \
 vode.o dgefa.o dgesl.o dgbfa.o dgbsl.o -lblas -o "$OUT"
sha256sum "$OUT"
