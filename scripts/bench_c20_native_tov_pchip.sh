#!/usr/bin/env bash
# C20 frontend runner: prebuild once, then isolate non-reentrant DVODE by process.
# The unchanged C17 source-locked campaign owns the 40 densities and 4440 shots.
set -euo pipefail
if [[ $# -ne 2 || ! "$1" =~ ^[12]$ || "$2" != /* || -e "$2" ]]; then
  echo "usage: $0 1|2 /abs/new-output-directory" >&2
  exit 2
fi
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$ROOT/python_frontend"
# Avoid the inherited cold-build shared temporary path in concurrent workers.
(cd "$ROOT"; python3 -c 'from compose_eos.stellar_cowling_native import _native_library; print("C20_PREBUILT",_native_library()._name)')
exec bash "$ROOT/scripts/bench_c17_extended_parallel.sh" "$1" "$2"
