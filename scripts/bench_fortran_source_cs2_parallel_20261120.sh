#!/usr/bin/env bash
# Independent Fortran-DVODE + source-EOS workers; legacy DVODE COMMON
# and SAVEd coefficient buffers are never shared among worker processes.
set -euo pipefail
if [[ $# -ne 2 ]];then
  echo "usage: $0 workers /abs/new-output" >&2; exit 2
fi
WORKERS=$1;OUT=$2
[[ "$WORKERS" =~ ^[12]$ ]] || exit 2
[[ "$OUT" = /* && ! -e "$OUT" ]] || exit 2
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
mkdir -p "$OUT"
export PYTHONPATH="$ROOT/python_frontend"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
start=$(date +%s%N)
printf '%s\n' M2 M3 | xargs -P "$WORKERS" -I '{}' \
 bash -c '/usr/bin/time -f "STAR_MAX_RSS_KB %M" \
   python3 "$1/scripts/run_fortran_source_cs2_case_20261120.py" \
   --case "$2" --output "$3/$2.json" >"$3/$2.log" 2>"$3/$2.memory"' \
 _ "$ROOT" '{}' "$OUT"
end=$(date +%s%N)
for model in M2 M3;do
  [[ -s "$OUT/$model.json" ]] || exit 3
  grep -q 'SOURCE_NATIVE_CASE' "$OUT/$model.log" || exit 3
done
(cd "$OUT";sha256sum M2.json M3.json > SCIENCE.sha256;
 sha256sum -c SCIENCE.sha256 >/dev/null)
awk -v p="$WORKERS" -v ns=$((end-start)) \
 'BEGIN {printf "REAL_SOURCE_NATIVE_WORKERS %d STARS 2 WALL_S %.6f\n",p,ns/1e9}' \
 |tee "$OUT/TIMING.txt"
(cd "$OUT";sha256sum SCIENCE.sha256)
