#!/usr/bin/env bash
# Source-locked independent stellar cases, each process owns native DVODE state.
set -euo pipefail
if [[ $# -ne 4 ]]; then
  echo "usage: $0 workers density-manifest output-directory frequencies-per-star" >&2
  exit 2
fi
WORKERS=$1; INPUT=$2; OUT=$3; COUNT=$4
[[ "$WORKERS" =~ ^[1-4]$ ]] || exit 2
[[ "$COUNT" =~ ^[1-9][0-9]*$ ]] || exit 2
(( COUNT >= 8 && COUNT <= 4096 )) || exit 2
[[ -f "$INPUT" ]] || exit 2
mapfile -t DENSITIES < "$INPUT"
(( ${#DENSITIES[@]} >= WORKERS && ${#DENSITIES[@]} <= 16 )) || exit 2
for d in "${DENSITIES[@]}"; do
  [[ "$d" =~ ^[0-9]+\.[0-9]+$ ]] || { echo "invalid density: $d" >&2; exit 2; }
done
[[ $(printf '%s\n' "${DENSITIES[@]}" | LC_ALL=C sort -u | wc -l) == "${#DENSITIES[@]}" ]] || {
  echo "density duplicates would overwrite a result" >&2; exit 2;
}
[[ ! -e "$OUT" ]] || { echo "output folder already exists" >&2; exit 2; }
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
mkdir -p "$OUT"
export PYTHONPATH="$ROOT/python_frontend"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
START=$(date +%s%N)
printf '%s\n' "${DENSITIES[@]}" | xargs -P "$WORKERS" -I '{}' \
  bash -c 'exec python3 "$1/scripts/run_cowling_native_bulk_scan_c9.py" \
    --central-density "$2" --output "$3/star_$2.json" --count "$4"' \
   _ "$ROOT" '{}' "$OUT" "$COUNT" > "$OUT/WORKERS.log"
END=$(date +%s%N)
for d in "${DENSITIES[@]}"; do
  [[ -s "$OUT/star_$d.json" ]] || { echo "missing result: $d" >&2; exit 3; }
  python3 -m json.tool "$OUT/star_$d.json" >/dev/null
done
(cd "$OUT"; sha256sum star_*.json > SCIENCE.sha256; sha256sum -c SCIENCE.sha256 >/dev/null)
awk -v ns=$((END-START)) -v n="$WORKERS" -v jobs="${#DENSITIES[@]}" \
  'BEGIN {printf "REAL_STAR_WORKERS %d CASES %d WALL_S %.6f\n",n,jobs,ns/1e9}' \
  | tee "$OUT/TIMING.txt"
(cd "$OUT"; sha256sum SCIENCE.sha256)
