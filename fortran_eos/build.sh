#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$ROOT/fortran_eos"
MODE=force
case "${1:-}" in
  --ensure) MODE=ensure; shift ;;
  --signature) MODE=signature; shift ;;
esac
BUILD="${1:-$ROOT/build/compose_eos_fortran}"
FC="${FC:-gfortran}"
FLAGS="${FFLAGS:--std=f2018 -O2 -ffree-line-length-none -fimplicit-none -fcheck=all}"
LEGACY_FLAGS="${SLSQP_FFLAGS:--std=legacy -O2 -ffixed-line-length-none}"
LINK_FLAGS="${LDFLAGS:-}"
NATIVE_LIBS="${NATIVE_LIBS:--llapack -lblas}"
COMPILER="$(command -v "$FC")"
signature() {
  {
    printf '%s\n' "compiler=$COMPILER" "FFLAGS=$FLAGS" "SLSQP_FFLAGS=$LEGACY_FLAGS" "LDFLAGS=$LINK_FLAGS" "NATIVE_LIBS=$NATIVE_LIBS"
    "$FC" --version
    "$FC" -dumpmachine
    sha256sum "$COMPILER" "$SRC/build.sh" "$SRC/run-generator.sh"
    find "$SRC/src" "$SRC/app" "$SRC/vendor" -type f \( -name '*.f90' -o -name '*.f' -o -name 'PROVENANCE.json' -o -name '*LICENSE*' \) -print0 |
      LC_ALL=C sort -z | xargs -0 sha256sum
  } | sha256sum | cut -d' ' -f1
}
if [[ "$MODE" == signature ]]; then
  signature
  exit 0
fi
mkdir -p "$BUILD"
BUILD="$(cd "$BUILD" && pwd)"
exec 9>"$BUILD/.build.lock"
flock 9
BEFORE="$(signature)"
if [[ "$MODE" == ensure && -x "$BUILD/compose_eos_fortran" && -x "$BUILD/compose_eos_generate" && -x "$BUILD/compose_eos_import_cold_table" &&
      -f "$BUILD/.source-signature" && "$(cat "$BUILD/.source-signature")" == "$BEFORE" ]]; then
  printf '%s\n' "$BUILD/compose_eos_fortran" "$BUILD/compose_eos_generate"
  exit 0
fi
STAGE="$(mktemp -d "$BUILD/.compile.XXXXXXXX")"
trap 'rm -rf "$STAGE"' EXIT
mkdir -p "$STAGE/mod" "$STAGE/obj"
read -r -a FARGS <<< "$FLAGS"
read -r -a LEGACY_ARGS <<< "$LEGACY_FLAGS"
read -r -a LINK_ARGS <<< "$LINK_FLAGS"
read -r -a LIB_ARGS <<< "$NATIVE_LIBS"
compile() {
  local name="$1" directory="${2:-src}"
  "$FC" "${FARGS[@]}" -J"$STAGE/mod" -I"$STAGE/mod" -c "$SRC/$directory/$name.f90" -o "$STAGE/obj/$name.o"
}
OLD_MODULES=(
  compose_eos_types compose_eos_rmf compose_eos_cld compose_eos_tov compose_eos_component
  compose_eos_css compose_eos_phase_transition compose_eos_matching compose_eos_reconstruct
  compose_eos_causality compose_eos_io compose_eos_verify compose_eos_wigner_seitz
)
OLD_OBJECTS=()
for name in "${OLD_MODULES[@]}"; do
  compile "$name"
  OLD_OBJECTS+=("$STAGE/obj/$name.o")
done
compile compose_eos_fortran app
# Original Netlib DVODE, separate-process legacy COMMON workspace (not thread safe).
"$FC" "${LEGACY_ARGS[@]}" -c "$SRC/vendor/odepack/vode.f" -o "$STAGE/obj/vode.o"
LINPACK_OBJECTS=()
for name in dgefa dgesl dgbfa dgbsl; do
  "$FC" "${LEGACY_ARGS[@]}" -c "$SRC/vendor/odepack/$name.f" -o "$STAGE/obj/$name.o"
  LINPACK_OBJECTS+=("$STAGE/obj/$name.o")
done
"$FC" "${FARGS[@]}" "${OLD_OBJECTS[@]}" "$STAGE/obj/compose_eos_fortran.o" "$STAGE/obj/vode.o" "${LINPACK_OBJECTS[@]}" "${LINK_ARGS[@]}" "${LIB_ARGS[@]}" -o "$STAGE/compose_eos_fortran"
"$FC" "${LEGACY_ARGS[@]}" -c "$SRC/vendor/slsqp/slsqp_optmz.f" -o "$STAGE/obj/slsqp_optmz.o"
NATIVE_MODULES=(
  compose_eos_optimize_native compose_eos_core_native compose_eos_core_domain compose_eos_gdfm_native
  compose_eos_lower_crust_native compose_eos_crust_attachment compose_eos_tf_models compose_eos_spg_native compose_eos_spatial_tf
  compose_eos_bps_primary compose_eos_bbp_primary compose_eos_nv_primary
)
NATIVE_OBJECTS=()
for name in "${NATIVE_MODULES[@]}"; do
  compile "$name"
  NATIVE_OBJECTS+=("$STAGE/obj/$name.o")
done
compile compose_eos_generate app
"$FC" "${FARGS[@]}" "$STAGE/obj/compose_eos_types.o" "$STAGE/obj/compose_eos_io.o" \
  "${NATIVE_OBJECTS[@]}" "$STAGE/obj/slsqp_optmz.o" "$STAGE/obj/compose_eos_generate.o" \
  "${LINK_ARGS[@]}" "${LIB_ARGS[@]}" -o "$STAGE/compose_eos_generate"
compile compose_eos_import_cold_table app
"$FC" "${FARGS[@]}" "$STAGE/obj/compose_eos_types.o" "$STAGE/obj/compose_eos_io.o" \
  "$STAGE/obj/compose_eos_import_cold_table.o" "${LINK_ARGS[@]}" -o "$STAGE/compose_eos_import_cold_table"
AFTER="$(signature)"
if [[ "$BEFORE" != "$AFTER" ]]; then
  printf '%s\n' 'Fortran sources or build configuration changed during compilation; refusing to publish a stale build.' >&2
  exit 75
fi
# Publish only after all executables have compiled successfully.
mkdir -p "$BUILD/mod" "$BUILD/obj"
cp "$STAGE/mod/"*.mod "$BUILD/mod/"
cp "$STAGE/obj/"*.o "$BUILD/obj/"
mv "$STAGE/compose_eos_fortran" "$BUILD/compose_eos_fortran"
mv "$STAGE/compose_eos_generate" "$BUILD/compose_eos_generate"
mv "$STAGE/compose_eos_import_cold_table" "$BUILD/compose_eos_import_cold_table"
printf '%s\n' "$AFTER" > "$BUILD/.source-signature.tmp"
mv "$BUILD/.source-signature.tmp" "$BUILD/.source-signature"
printf '%s\n' "$BUILD/compose_eos_fortran" "$BUILD/compose_eos_generate"
