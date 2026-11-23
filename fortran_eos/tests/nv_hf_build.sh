#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD="${1:-$ROOT/build/nv_hf_tests}"
mkdir -p "$BUILD/mod"
FC="${FC:-gfortran}"
FFLAGS="${FFLAGS:--std=f2018 -O2 -ffree-line-length-none -fimplicit-none -fcheck=all -ffpe-trap=invalid,zero,overflow}"
read -r -a FLAGS <<< "$FFLAGS"
SOURCES=(compose_eos_types compose_eos_nv_primary compose_eos_nv_original compose_eos_nv_dme compose_eos_nv_hf)
OBJECTS=()
for name in "${SOURCES[@]}"; do
 "$FC" "${FLAGS[@]}" -J"$BUILD/mod" -I"$BUILD/mod" -c "$ROOT/fortran_eos/src/$name.f90" -o "$BUILD/$name.o"
 OBJECTS+=("$BUILD/$name.o")
done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for name in nv_hf_numerics nv_hf_variation nv_hf_basis_extension; do
 "$FC" "${FLAGS[@]}" -J"$BUILD/mod" -I"$BUILD/mod" "${OBJECTS[@]}" "$ROOT/fortran_eos/tests/$name.f90" -llapack -lblas -o "$BUILD/$name"
 "$BUILD/$name" "$ROOT/data/eos_original_papers_20261123/research_nv/nv1973_original_source.dat" | tee "$BUILD/$name.log"
done
sha256sum "$BUILD/nv_hf_numerics" "$BUILD/nv_hf_variation" "$BUILD/nv_hf_basis_extension" > "$BUILD/SHA256SUMS.txt"
