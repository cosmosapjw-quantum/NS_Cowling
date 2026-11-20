#!/usr/bin/env bash
# Reproducible Fortran-native benchmark, including a parent-source control.
set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "usage: $0 /abs/path/to/native-probe current|parent" >&2
  exit 2
fi
OUT=$1
MODE=$2
[[ "$OUT" = /* ]] || { echo "output path must be absolute" >&2; exit 2; }
[[ "$MODE" == current || "$MODE" == parent ]] || exit 2
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
DIR=$(dirname "$OUT")
mkdir -p "$DIR"
if [[ "$MODE" == parent ]]; then
  CORE="$DIR/compose_eos_cowling_dvode.parent.f90"
  git -C "$ROOT" show \
    d6fdd55d271f397dc7fc9164cddc416ed9ecf655:fortran_eos/src/compose_eos_cowling_dvode.f90 \
    > "$CORE"
else
  CORE="$ROOT/fortran_eos/src/compose_eos_cowling_dvode.f90"
fi
gfortran -O2 -ffree-line-length-none -J "$DIR" \
  "$CORE" \
  "$ROOT/fortran_eos/vendor/odepack/vode.f" \
  "$ROOT/fortran_eos/vendor/odepack/dgefa.f" \
  "$ROOT/fortran_eos/vendor/odepack/dgesl.f" \
  "$ROOT/fortran_eos/vendor/odepack/dgbfa.f" \
  "$ROOT/fortran_eos/vendor/odepack/dgbsl.f" \
  "$ROOT/fortran_eos/tests/cowling_c8_native_probe.f90" \
  -lblas -o "$OUT" 2>"$OUT.build.log"
echo "BUILT $MODE $OUT"
sha256sum "$CORE" "$OUT"
