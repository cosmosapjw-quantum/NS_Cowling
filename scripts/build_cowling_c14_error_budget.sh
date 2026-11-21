#!/usr/bin/env bash
# C14 Fortran-only original-EOS eigenfrequency 3-axis error-budget executable.
# Build project sources with debug bounds or optimized flags; frozen Netlib
# F77 compatibility remains separately compiled, not certified bounds-safe.
set -euo pipefail
if [[ $# -ne 2 ]];then
  echo "usage: $0 /abs/output-executable checked|optimized" >&2
  exit 2
fi
OUT=$1;MODE=$2
[[ "$OUT" = /* ]] || exit 2
[[ "$MODE" == checked || "$MODE" == optimized ]] || exit 2
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
DIR=$(dirname "$OUT")
mkdir -p "$DIR"
cd "$DIR"
FLAGS=(-O0 -g -fcheck=all -ffree-line-length-none)
if [[ "$MODE" == optimized ]];then
  FLAGS=(-O2 -ffree-line-length-none)
fi
gfortran "${FLAGS[@]}" -c \
 "$ROOT/fortran_eos/src/compose_eos_cowling_dvode.f90" \
 "$ROOT/fortran_eos/src/compose_eos_cowling_pchip.f90" \
 "$ROOT/fortran_eos/src/compose_eos_cowling_source_exact.f90" \
 "$ROOT/fortran_eos/src/compose_eos_cowling_sharp_source_root.f90" \
 "$ROOT/fortran_eos/src/compose_eos_cowling_error_budget.f90" \
 "$ROOT/fortran_eos/tests/cowling_c14_error_budget.f90"
gfortran -O2 -std=legacy -c \
 "$ROOT/fortran_eos/vendor/odepack/vode.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgefa.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgesl.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgbfa.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgbsl.f"
gfortran *.o -lblas -o "$OUT"
sha256sum "$OUT"
