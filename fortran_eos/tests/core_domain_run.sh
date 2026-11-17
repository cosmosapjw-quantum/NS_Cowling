#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="${1:-$ROOT/data/research_blockers_20261116/core}"
mkdir -p "$OUT/build"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
gfortran -std=f2018 -O2 -Wall -Wextra -ffree-line-length-none -fcheck=all \
  -J"$OUT/build" -I"$OUT/build" \
  "$ROOT/fortran_eos/src/compose_eos_types.f90" \
  "$ROOT/fortran_eos/src/compose_eos_core_native.f90" \
  "$ROOT/fortran_eos/src/compose_eos_core_domain.f90" \
  "$ROOT/fortran_eos/tests/core_domain_probe.f90" \
  -o "$OUT/build/core_domain_probe" 2>"$OUT/compile.log"
cd "$OUT"
set +e
./build/core_domain_probe >run.log 2>&1
rc=$?
set -e
printf 'driver_exit_code=%d\n' "$rc" >exit_status.txt
exit "$rc"
