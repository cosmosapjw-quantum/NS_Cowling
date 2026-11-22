#!/usr/bin/env bash
# C17 extended source-locked original EOS f/p1 density sequence with overlap gate.
# Fortran DVODE COMMON is process-local and never entered by concurrent threads.
set -euo pipefail
if [[ $# -ne 2 ]];then
  echo "usage: $0 1|2 /abs/new-output-directory" >&2;exit 2
fi
WORKERS=$1;OUT=$2
[[ "$WORKERS" =~ ^[12]$ ]] || exit 2
[[ "$OUT" = /* && ! -e "$OUT" ]] || exit 2
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
mkdir -p "$OUT"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$ROOT/python_frontend"
start=$(date +%s%N)
printf '%s\n' M2 M3 | xargs -P "$WORKERS" -I '{}' \
 bash -c '/usr/bin/time -f "C17_PEAK_RSS_KB %M" \
 python3 "$1/scripts/run_c17_extended_density.py" \
 --case "$2" --output "$3/$2.json" \
 >"$3/$2.log" 2>"$3/$2.memory"' \
 _ "$ROOT" '{}' "$OUT"
end=$(date +%s%N)
for model in M2 M3;do
 [[ -s "$OUT/$model.json" ]] || exit 3
 grep -q "C17_WORKER_DONE $model" "$OUT/$model.log"
 [[ -s "$OUT/${model}_all_scan_points.json.gz" ]] || exit 3
 python3 -m json.tool "$OUT/$model.json" >/dev/null
done
(cd "$OUT";sha256sum M2.json M3.json \
    M2_all_scan_points.json.gz M3_all_scan_points.json.gz >SCIENTIFIC.sha256;
 sha256sum -c SCIENTIFIC.sha256 >/dev/null)
awk -v p="$WORKERS" -v n=$((end-start)) \
 'BEGIN {printf "C17_PROCESSES %d STARS 2 DENSITIES 40 WALL_S %.6f\n",p,n/1e9}' |tee "$OUT/TIMING.txt"
(cd "$OUT";sha256sum SCIENTIFIC.sha256)
