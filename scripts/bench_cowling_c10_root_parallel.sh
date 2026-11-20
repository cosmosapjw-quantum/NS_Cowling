#!/usr/bin/env bash
# Run one independent Fortran TOV + Fortran Cowling scan/root per process.
# No OpenMP sharing: legacy Netlib DVODE has process-wide COMMON state.
set -euo pipefail
if [[ $# -ne 3 ]]; then
  echo "usage: $0 workers central-density-file new-output-dir" >&2;exit 2
fi
WORKERS=$1;INPUT=$2;OUT=$3
[[ "$WORKERS" =~ ^[1-4]$ ]] || exit 2
[[ -f "$INPUT" && ! -e "$OUT" ]] || exit 2
mapfile -t DS <"$INPUT"
(( ${#DS[@]}>=WORKERS && ${#DS[@]}<=16 )) || exit 2
for d in "${DS[@]}";do
  [[ "$d" =~ ^[0-9]+\.[0-9]+$ ]] || exit 2
done
[[ $(printf '%s\n' "${DS[@]}" | sort -u | wc -l) == "${#DS[@]}" ]] || exit 2
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
mkdir -p "$OUT"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
start=$(date +%s%N)
printf '%s\n' "${DS[@]}" | xargs -P "$WORKERS" -I '{}' \
  bash -c 'exec python3 "$1/scripts/run_cowling_c10_native_root_case.py" \
   --central-density "$2" --output "$3/star_$2.json"' \
   _ "$ROOT" '{}' "$OUT" > "$OUT/WORKERS.log"
stop=$(date +%s%N)
for d in "${DS[@]}";do
  [[ -s "$OUT/star_$d.json" ]] || exit 3
  python3 -m json.tool "$OUT/star_$d.json" >/dev/null
done
(cd "$OUT"; sha256sum star_*.json >SCIENTIFIC.sha256;sha256sum -c SCIENTIFIC.sha256 >/dev/null)
awk -v ns=$((stop-start)) -v workers="$WORKERS" -v jobs="${#DS[@]}" \
 'BEGIN {printf "C10_ROOT_PROCESSES %d STARS %d WALL_SECONDS %.6f\n",workers,jobs,ns/1e9}' \
 | tee "$OUT/TIMING.txt"
cat "$OUT/WORKERS.log"
(cd "$OUT";sha256sum SCIENTIFIC.sha256)
