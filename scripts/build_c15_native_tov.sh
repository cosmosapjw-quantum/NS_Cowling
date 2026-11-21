#!/usr/bin/env bash
# C15 independent two-phase native Fortran TOV source-locked replay builder.
# Project Fortran uses bounded compilation; untouched Netlib ODEPACK uses
# legacy-compatible -O2, which does NOT certify vendor-wide bounds safety.
set -euo pipefail
if [[ $# -ne 2 ]];then
 echo "usage: $0 /abs/new-executable checked|optimized" >&2;exit 2
fi
OUT=$1;MODE=$2
[[ "$OUT" = /* ]] || exit 2
[[ "$MODE" == checked || "$MODE" == optimized ]] || exit 2
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
DIR=$(dirname "$OUT")
mkdir -p "$DIR"
cd "$DIR"
FLAGS=(-O0 -g -fcheck=all -ffree-line-length-none)
if [[ "$MODE" == optimized ]];then FLAGS=(-O2 -ffree-line-length-none);fi
gfortran "${FLAGS[@]}" -c \
 "$ROOT/fortran_eos/src/compose_eos_cowling_pchip.f90" \
 "$ROOT/fortran_eos/src/compose_eos_cowling_tov_dvode.f90" \
 "$ROOT/fortran_eos/src/compose_eos_cowling_tov_two_phase_sweep.f90" \
 "$ROOT/fortran_eos/tests/cowling_c15_original_two_phase_tov.f90"
gfortran -O2 -std=legacy -c \
 "$ROOT/fortran_eos/vendor/odepack/vode.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgefa.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgesl.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgbfa.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgbsl.f"
gfortran *.o -lblas -o "$OUT"
sha256sum "$OUT"
