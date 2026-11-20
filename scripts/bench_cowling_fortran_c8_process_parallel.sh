#!/usr/bin/env bash
# Safe coarse-grained parallelism: independent processes only.
# Netlib DVODE uses COMMON and the Fortran modules keep SAVEd tables.
# No concurrent calls or OpenMP threads enter an individual DVODE image.
set -euo pipefail
if [[ $# -ne 5 ]]; then
  echo "usage: $0 /abs/native_probe workers jobs shots /abs/output-dir" >&2
  exit 2
fi
BINARY=$1
WORKERS=$2
JOBS=$3
SHOTS=$4
OUT=$5
[[ -x "$BINARY" ]] || { echo "missing executable" >&2; exit 2; }
[[ "$WORKERS" =~ ^[1-9][0-9]*$ && "$JOBS" =~ ^[1-9][0-9]*$ && "$SHOTS" =~ ^[1-9][0-9]*$ ]] || exit 2
(( WORKERS <= JOBS && JOBS <= 32 && SHOTS <= 100000 )) || exit 2
mkdir -p "$OUT"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
# Benchmark timing includes process launch, output serialization and teardown.
start_ns=$(date +%s%N)
seq 1 "$JOBS" | xargs -P "$WORKERS" -I '{}' sh -c \
  'exec "$1" "$2" "$3" > "$4/job_$2.log"' _ "$BINARY" '{}' "$SHOTS" "$OUT"
end_ns=$(date +%s%N)
for ((j=1;j<=JOBS;j++)); do
  [[ $(grep -c '^C8_JOB' "$OUT/job_$j.log") == 1 ]] || {
    echo "missing native result for $j" >&2; exit 3;
  }
done
awk '$1 == "C8_JOB" {print}' "$OUT"/job_*.log | LC_ALL=C sort -k2,2n > "$OUT/SCIENTIFIC_RESULTS.txt"
sha256sum "$OUT/SCIENTIFIC_RESULTS.txt" > "$OUT/SCIENTIFIC_RESULTS.sha256"
awk -v elapsed=$((end_ns-start_ns)) -v workers="$WORKERS" \
  -v jobs="$JOBS" -v shots="$SHOTS" \
  'BEGIN { printf "PARALLEL_WORKERS %d JOBS %d SHOTS_PER_JOB %d WALL_S %.6f\n", workers,jobs,shots,elapsed/1e9 }' | tee "$OUT/TIMING.txt"
cat "$OUT/SCIENTIFIC_RESULTS.sha256"
