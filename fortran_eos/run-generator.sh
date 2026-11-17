#!/usr/bin/env bash
set -euo pipefail
if [[ "$#" != 2 ]]; then
  printf '%s\n' 'usage: run-generator.sh REQUEST_NML OUTPUT_DIRECTORY' >&2
  exit 2
fi
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="${COMPOSE_EOS_NATIVE_BUILD_DIR:-$ROOT/build/compose_eos_fortran}"
bash "$ROOT/fortran_eos/build.sh" --ensure "$BUILD" >&2
exec "$BUILD/compose_eos_generate" "$1" "$2"
