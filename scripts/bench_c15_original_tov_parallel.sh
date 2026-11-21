#!/usr/bin/env bash
# C15 original EOS, TOV background and Cowling p1 source-locked process sweep.
# One independent OS worker per stellar model, DVODE COMMON is not thread-safe.
set -euo pipefail
if [[ $# -ne 2 ]];then
  echo "usage: $0 1|2 /abs/new-output-directory" >&2
  exit 2
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
 bash -c '/usr/bin/time -f "MAX_RSS_KB %M" \
  python3 "$1/scripts/bench_c15_original_tov_uncertainty.py" \
  --case "$2" --output "$3/$2.json" --quiet \
  >"$3/$2.log" 2>"$3/$2.memory"' \
 _ "$ROOT" '{}' "$OUT"
end=$(date +%s%N)
for key in M2 M3;do
 [[ -s "$OUT/$key.json" ]] || exit 3
 grep -q "^C15_FINAL $key" "$OUT/$key.log"
 python3 -m json.tool "$OUT/$key.json" >/dev/null
done
(cd "$OUT";sha256sum M2.json M3.json >SCIENTIFIC.sha256;
 sha256sum -c SCIENTIFIC.sha256 >/dev/null)
awk -v workers="$WORKERS" -v secs=$((end-start)) \
 'BEGIN{printf "C15_WORKERS %d TOV_CONTROLS 40 WALL_S %.6f\n",workers,secs/1e9}' |tee "$OUT/TIMING.txt"
(cd "$OUT";sha256sum SCIENTIFIC.sha256)
