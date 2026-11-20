#!/usr/bin/env bash
# Original source-locked two-star pure Fortran process-parallel regression.
# A worker repeats Adams and BDF roots N times and rejects any byte-drift.
# No thread enters DVODE's global COMMON or shared coefficient arrays.
set -euo pipefail
if [[ $# -ne 5 ]];then
  echo "usage: $0 /abs/fortran-binary /abs/input-dir workers repeats /abs/new-output" >&2
  exit 2
fi
BIN=$1; INPUT=$2;WORKERS=$3;REPEATS=$4;OUT=$5
[[ -x "$BIN" && -d "$INPUT" ]] || exit 2
[[ "$WORKERS" =~ ^[12]$ ]] || exit 2
[[ "$REPEATS" =~ ^[1-9][0-9]*$ ]] || exit 2
(( REPEATS<=1000 )) || exit 2
[[ "$OUT" = /* && ! -e "$OUT" ]] || exit 2
for model in M2 M3;do
  [[ -s "$INPUT/ORIGINAL_FORTRAN_ROOT_$model.dat" ]] || exit 2
done
mkdir -p "$OUT"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
start=$(date +%s%N)
printf '%s\n' M2 M3 | xargs -P "$WORKERS" -I '{}' \
  bash -c '/usr/bin/time -f "C12_MAX_RSS_KB %M" \
    "$1" "$2/ORIGINAL_FORTRAN_ROOT_$3.dat" "$4" \
    >"$5/$3.log" 2>"$5/$3.memory"' \
  _ "$BIN" "$INPUT" '{}' "$REPEATS" "$OUT"
end=$(date +%s%N)
for model in M2 M3;do
  [[ -s "$OUT/$model.log" ]] || exit 3
  grep -q 'C12_FORTRAN_SOURCE_EXACT_SHARP_ROOT_ADAMS_BDF_PASS' "$OUT/$model.log"
  grep -q "C12_ROOT_REPEAT_DETERMINISTIC *$REPEATS" "$OUT/$model.log"
  [[ $(grep -c '^C12_NATIVE_SHARP_ROOT' "$OUT/$model.log") == 2 ]] || exit 3
done
(cd "$OUT";sha256sum M2.log M3.log > SCIENCE.sha256;
  sha256sum -c SCIENCE.sha256 >/dev/null)
awk -v n="$WORKERS" -v ns=$((end-start)) -v rep="$REPEATS" \
 'BEGIN {printf "C12_FORTRAN_WORKERS %d ORIGINAL_STARS 2 ROOT_SHOTS_PER_STAR %d WALL_S %.6f\n", n, 2*rep, ns/1e9}' \
 | tee "$OUT/TIMING.txt"
(cd "$OUT";sha256sum SCIENCE.sha256)
