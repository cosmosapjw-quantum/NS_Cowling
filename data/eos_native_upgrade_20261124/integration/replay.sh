#!/usr/bin/env bash
# All physical solvers under test run in Fortran. Use a fresh scratch directory.
set -euo pipefail
if [[ $# != 4 ]]; then
  echo "usage: $0 FULL_BUILD GENERATED_CATALOG_DIR ARCHIVED_OR_PINNED_RUNS EMPTY_SCRATCH" >&2
  exit 2
fi
repo=$(cd "$(dirname "$0")/../../.." && pwd)
full=$(realpath "$1")
catalog=$(realpath "$2")
baseline=$(realpath "$3")
mkdir -p "$4"
out=$(realpath "$4")
[[ -z $(ls -A "$out") ]] || { echo "scratch directory must be empty" >&2; exit 2; }
umask 077
mkdir -p "$out"/{obj,checks,private,tmp}
exec 9<"$full/.build.lock"
flock -s 9
export COMPOSE_EOS_FORTRAN_EXE="$full/compose_eos_fortran"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$repo/python_frontend${PYTHONPATH:+:$PYTHONPATH}"
fc=${FC:-gfortran}
read -r -a libs <<< "${NATIVE_LIBS:--llapack -lblas}"
flags=(-std=f2018 -O2 -ffree-line-length-none -fimplicit-none -fcheck=all -fopenmp)
tests="$repo/fortran_eos/tests"
obj="$full/obj"
vendor=("$obj"/{vode,dgefa,dgesl,dgbfa,dgbsl}.o)
cd "$out"
sha256sum "$full/compose_eos_fortran" | awk '{print $1}' > binary.sha256
cat "$full/.source-signature" > source-signature
"$fc" --version | head -1
printf 'full native executable %s\n' "$(cat binary.sha256)"
check() {
  local name=$1
  shift
  "$fc" "${flags[@]}" -I "$full/mod" -I "$out/obj" -J "$out/obj" \
    "$tests/$name.f90" "$@" "${libs[@]}" -o "$out/checks/$name"
}
run_check() {
  local name=$1
  shift
  "$out/checks/$name" "$@" > "$out/checks/$name.log" 2>&1 || { cat "$out/checks/$name.log"; return 1; }
  printf 'PASS %s\n' "$name"
}
for name in test_matching_{slice_conventions,search_exhaustive_20261123}; do
  check "$name" "$obj/compose_eos_matching.o" "$obj/compose_eos_types.o"
  run_check "$name"
done
for name in test_tov_love test_tov_love_literature_20261118 test_tidal_rhs_literature_20261118 test_tov_layered; do
  check "$name" "$obj"/compose_eos_{tov,verify,types}.o "${vendor[@]}"
  run_check "$name"
done
for module in dvode sharp_root; do
  "$fc" "${flags[@]}" -I "$out/obj" -J "$out/obj" -c \
    "$repo/fortran_eos/src/compose_eos_cowling_$module.f90" -o "$out/obj/$module.o"
done
for name in cowling_c9_bulk_scan_test cowling_c10_brent_root_test; do
  check "$name" "$out/obj/dvode.o" "${vendor[@]}"
  run_check "$name"
done
ref=cab027e9b5ccbc42d6509531f7a6201f2c13d184
sharp=data/cowling_sharp_native_root_20261120/ORIGINAL_SOURCE_PCHIP_ONEJUMP.dat.gz
pk1=data/native_dense_interpolation_20261118/rmf/PK1_CANONICAL/table/compose
git -C "$repo" archive "$ref" "$sharp" "$pk1" | tar xf - -C "$out"
gzip -dc "$out/$sharp" > "$out/onejump.dat"
check cowling_sharp_onejump_root_test "$out/obj/sharp_root.o" "$out/obj/dvode.o" "${vendor[@]}"
run_check cowling_sharp_onejump_root_test "$out/onejump.dat"
# The existing CLI fixture resolves its EOS directory relative to its own location.
mkdir -p "$out/python_frontend/tests"
cp "$repo/python_frontend/tests/test_native_tov_batch_20261118.py" "$out/python_frontend/tests/"
python -m pytest -q -p no:cacheprovider --noconftest "$out/python_frontend/tests/test_native_tov_batch_20261118.py"
python "$repo/scripts/run_eos_all_pairs_20261123.py" \
  --catalog "$catalog/legacy_control.json" --native-exe "$full/compose_eos_fortran" \
  --out "$out/runs/legacy_control" --private-artifacts "$out/private/legacy_control" \
  --work "$out/tmp" --n-bridge 80 --workers 1 > "$out/legacy.log"
python "$tests/compare_matching_replays.py" "$baseline" "$out/runs" "$out/legacy_comparison.json" legacy_control
# The archived harness explicitly exposes private kernels in a test-local module copy.
python "$tests/run_matching_native_adversarial.py" --build-dir "$out/probe" --log-dir "$out/probe_logs" > "$out/probe.log"
printf 'PASS four native private-kernel controls\n'
