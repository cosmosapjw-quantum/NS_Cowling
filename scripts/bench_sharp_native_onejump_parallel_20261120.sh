#!/usr/bin/env bash
# Fortran-only numerical workers on one source-locked one-Maxwell PCHIP input.
# Shared COMMON DVODE and SAVEd arrays are kept in separate OS processes.
set -euo pipefail
if [[ $# -ne 4 ]]; then
  echo "usage: $0 /abs/native-driver /abs/input.dat workers /abs/new-output" >&2
  exit 2
fi
BIN=$1;INPUT=$2;WORKERS=$3;OUT=$4
[[ -x "$BIN" && -s "$INPUT" ]] || exit 2
[[ "$WORKERS" =~ ^[124]$ ]] || exit 2
[[ "$OUT" = /* && ! -e "$OUT" ]] || exit 2
mkdir -p "$OUT"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
start=$(date +%s%N)
seq 1 4 | xargs -P "$WORKERS" -I '{}' sh -c \
  '/usr/bin/time -f "MAX_RSS_KB %M" "$1" "$2" "$3" \
    > "$4/job_$3.log" 2> "$4/job_$3.memory"' \
    _ "$BIN" "$INPUT" '{}' "$OUT"
end=$(date +%s%N)
for j in 1 2 3 4; do
  [[ $(grep -c '^FORTRAN_SHARP_CASE' "$OUT/job_$j.log") == 2 ]] || {
    echo "worker $j produced an incomplete native result" >&2; exit 3;
  }
  grep -q 'SHARP_ONEJUMP_NATIVE_SELECTED_CASE_PASS' "$OUT/job_$j.log" || exit 3
done
awk '$1 == "FORTRAN_SHARP_CASE" {print}' "$OUT"/job_*.log |
  LC_ALL=C sort -k2,2n -k6,6n > "$OUT/SCIENTIFIC_RESULTS.txt"
(cd "$OUT";sha256sum SCIENTIFIC_RESULTS.txt > SCIENCE.sha256;
 sha256sum -c SCIENCE.sha256 >/dev/null)
awk -v w="$WORKERS" -v t=$((end-start)) \
 'BEGIN {printf "ONEJUMP_FORTRAN_WORKERS %d MODES 4 WALL_SECONDS %.6f\n",w,t/1e9}' \
 | tee "$OUT/TIMING.txt"
(cd "$OUT";sha256sum SCIENTIFIC_RESULTS.txt)
