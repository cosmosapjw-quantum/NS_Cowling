#!/usr/bin/env bash
# Run Fortran original-source C14 error-budget in isolated OS processes.
# The source EOS and geometry materialization are already frozen in C12.
# Timings start AFTER gzip input decompression.
set -euo pipefail
if [[ $# -ne 4 ]];then
 echo "usage: $0 /abs/fortran-exe /abs/source-input-dir 1|2 /abs/NEW-output" >&2;exit 2
fi
BIN=$1;INPUT=$2;WORKERS=$3;OUT=$4
[[ -x "$BIN" && -d "$INPUT" ]] || exit 2
[[ "$WORKERS" =~ ^[12]$ ]] || exit 2
[[ "$OUT" = /* && ! -e "$OUT" ]] || exit 2
mkdir -p "$OUT"
for key in M2 M3;do
  gzip -dc "$INPUT/ORIGINAL_FORTRAN_ROOT_$key.dat.gz" >"$OUT/$key.dat"
done
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
start=$(date +%s%N)
printf '%s\n' M2 M3 | xargs -P "$WORKERS" -I '{}' \
 bash -c '/usr/bin/time -f "C14_PEAK_RSS_KB %M" "$1" "$2/$3.dat" \
     >"$2/$3.science" 2>"$2/$3.memory"' \
 _ "$BIN" "$OUT" '{}'
end=$(date +%s%N)
for key in M2 M3;do
 grep -q 'C14_FORTRAN_3AXIS_ERROR_BUDGET_COMPLETE' "$OUT/$key.science"
 [[ $(grep -c 'C14_ROW' "$OUT/$key.science") == 24 ]] || {
 echo "C14 failed number of rows for $key" >&2;exit 3;}
done
(cd "$OUT";sha256sum M2.science M3.science >SCIENTIFIC.sha256;sha256sum -c SCIENTIFIC.sha256 >/dev/null)
# No duplicate uncompressed multi-MB EOS inputs remain in the evidence folder:
# each one is regenerated from the immutable, tracked C12 .dat.gz files.
for key in M2 M3;do rm -f -- "$OUT/$key.dat";done
awk -v p="$WORKERS" -v ns=$((end-start)) \
  'BEGIN {printf "C14_WORKERS %d ROOTS 48 WALL_S %.6f\n",p,ns/1e9}' |tee "$OUT/TIMING.txt"
(cd "$OUT";sha256sum SCIENTIFIC.sha256)
