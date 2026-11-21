#!/usr/bin/env bash
# C13 one physical original EOS star per independent OS process.
# DVODE COMMON arrays are NONREENTRANT, not exposed to shared threading.
set -euo pipefail
if [[ $# -ne 3 ]];then
  echo "usage: $0 worker_count samples /abs/new-output" >&2
  exit 2
fi
WORKERS=$1;SAMPLES=$2;OUT=$3
[[ "$WORKERS" =~ ^[12]$ ]] || exit 2
[[ "$SAMPLES" =~ ^[1-9][0-9]*$ ]] || exit 2
[[ "$OUT" = /* && ! -e "$OUT" ]] || exit 2
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
mkdir -p "$OUT"
export PYTHONPATH="$ROOT/python_frontend"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
start=$(date +%s%N)
printf '%s\n' M2 M3 | xargs -P "$WORKERS" -I '{}' \
 bash -c '/usr/bin/time -f "MAX_RSS_KB %M" python3 \
   "$1/scripts/run_c13_native_pchip_star.py" --case "$2" \
   --samples "$3" --output "$4/$2.json" \
   > "$4/$2.log" 2>"$4/$2.memory"' \
 _ "$ROOT" '{}' "$SAMPLES" "$OUT"
end=$(date +%s%N)
for m in M2 M3;do
 [[ -s "$OUT/$m.json" ]] || exit 3
 grep -q '^C13_FORTRAN_STAR' "$OUT/$m.log" || exit 3
 python3 -m json.tool "$OUT/$m.json" >/dev/null
done
(cd "$OUT";sha256sum M2.json M3.json >SCIENTIFIC.sha256;sha256sum -c SCIENTIFIC.sha256 >/dev/null)
awk -v p="$WORKERS" -v ns=$((end-start)) \
 'BEGIN {printf "C13_NATIVE_WORKERS %d STARS 2 WALL_SECONDS %.6f\n",p,ns/1e9}' \
 | tee "$OUT/TIMING.txt"
(cd "$OUT";sha256sum SCIENTIFIC.sha256)
