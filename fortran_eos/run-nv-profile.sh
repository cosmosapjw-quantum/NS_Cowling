#!/usr/bin/env bash
# Usage: run-nv-profile.sh OUTPUT_DIRECTORY [RADIAL_CELLS] [DENSITY_NODES]
set -eu
root=$(cd "$(dirname "$0")/.." && pwd)
out=${1:?supply an output directory}
mkdir -p "$out"
out=$(cd "$out" && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
cd "$tmp"
mkdir -p tables/{printed,d03x100,uniform_beta}
sources=()
for name in types io nv_primary nv_original nv_dme nv_hf nv_uniform_branch; do
    sources+=("$root/fortran_eos/src/compose_eos_${name}.f90")
done
# LAPACK_LIBS may name versioned libraries where development symlinks are absent.
read -r -a lapack_options <<< "${LAPACK_LIBS:--llapack -lblas}"
"${FC:-gfortran}" -std=f2008 -ffree-line-length-none -O3 "${sources[@]}" \
   "$root/fortran_eos/app/nv_profile_generate.f90" "${lapack_options[@]}" -o nv-profile
./nv-profile "$root/data/eos_original_papers_20261123/research_nv/nv1973_original_source.dat" "$tmp/tables" "${2:-128}" "${3:-49}"
cp -R tables/. "$out/"
