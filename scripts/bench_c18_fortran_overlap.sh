#!/usr/bin/env bash
# Native Fortran-only, pinned-base benchmark of the phasewise W overlap.
# No DVODE parallelization: legacy VODE COMMON remains non-reentrant.
set -euo pipefail
if [[ $# -ne 1 || "${1:-}" != /* || -e "${1:-}" ]];then
 echo "usage: $0 /absolute/new-output-directory" >&2;exit 2
fi
OUT=$1
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
BASE=b956e0a21cc31623d5034ab9f0cf9153b46a1d41
mkdir -p "$OUT"/{before,after}
git -C "$ROOT" show "$BASE:fortran_eos/src/compose_eos_cowling_mode_overlap.f90" >"$OUT/before/source.f90"
cp "$ROOT/fortran_eos/src/compose_eos_cowling_mode_overlap.f90" "$OUT/after/source.f90"
for side in before after;do
 (
  cd "$OUT/$side"
  gfortran -O2 -fopenmp -ffree-line-length-none -c source.f90 \
    "$ROOT/fortran_eos/tests/cowling_c18_overlap_benchmark.f90"
  gfortran -O2 -fopenmp *.o -o benchmark
 )
done
echo "BASE=$BASE" |tee "$OUT/perf.log"
echo "GFORT=$(gfortran -dumpfullversion)" |tee -a "$OUT/perf.log"
sha256sum "$OUT"/{before,after}/source.f90 |tee -a "$OUT/perf.log"
for workers in 1 4;do
 for rep in 1 2 3 4 5;do
  if (( rep % 2 == 0 ));then order=(after before);else order=(before after);fi
  for side in "${order[@]}";do
   echo "SIDE=$side REP=$rep WORKERS=$workers" |tee -a "$OUT/perf.log"
   OMP_NUM_THREADS=$workers OPENBLAS_NUM_THREADS=1 \
      "$OUT/$side/benchmark" 20000 |tee -a "$OUT/perf.log"
  done
 done
done
echo "C18_MICROBENCH_COMPLETE" |tee -a "$OUT/perf.log"
