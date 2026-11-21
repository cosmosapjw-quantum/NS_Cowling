#!/usr/bin/env bash
# C16 pure-Fortran f/p1 source-exact scanner and native nodal/shape regression.
# Frozen Netlib ODEPACK F77 separately compiles with legacy flags.
set -euo pipefail
if [[ $# -ne 2 ]];then
 echo "usage: $0 /abs/new-output-executable checked|optimized" >&2;exit 2
fi
BIN=$1;MODE=$2
[[ "$BIN" = /* ]] || exit 2
[[ "$MODE" = checked || "$MODE" = optimized ]] || exit 2
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
mkdir -p "$(dirname "$BIN")"
cd "$(dirname "$BIN")"
flags=(-O0 -g -fcheck=all -ffree-line-length-none)
if [[ "$MODE" = optimized ]];then flags=(-O2 -ffree-line-length-none);fi
gfortran "${flags[@]}" -c \
 "$ROOT/fortran_eos/src/compose_eos_cowling_dvode.f90" \
 "$ROOT/fortran_eos/src/compose_eos_cowling_source_exact.f90" \
 "$ROOT/fortran_eos/src/compose_eos_cowling_sharp_source_root.f90" \
 "$ROOT/fortran_eos/src/compose_eos_cowling_source_scan.f90" \
 "$ROOT/fortran_eos/src/compose_eos_cowling_mode_nodes.f90" \
 "$ROOT/fortran_eos/src/compose_eos_cowling_mode_overlap.f90" \
 "$ROOT/fortran_eos/tests/cowling_c16_original_source_f_p1.f90"
gfortran -O2 -std=legacy -c \
 "$ROOT/fortran_eos/vendor/odepack/vode.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgefa.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgesl.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgbfa.f" \
 "$ROOT/fortran_eos/vendor/odepack/dgbsl.f"
gfortran *.o -lblas -o "$BIN"
sha256sum "$BIN"
