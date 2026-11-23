#!/usr/bin/env bash
# Usage: OMP_NUM_THREADS=4 run-nv-profile.sh OUTPUT_DIRECTORY [RADIAL_CELLS] [DENSITY_NODES]
# Parallelize candidate evaluation inside each density, preserving density continuation.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
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
compile_flags=(-std=f2008 -ffree-line-length-none -O3 -fopenmp)
link_flags=()
default_libraries="-llapack -lblas"
if [[ $(uname -s) == Linux ]]; then
    # Only the orbital solver in the shared HF module uses LAPACK. Keep the
    # radial field routines used by this executable and drop the unused solver.
    compile_flags+=(-ffunction-sections -fdata-sections)
    link_flags+=(-Wl,--gc-sections)
    default_libraries=""
fi
read -r -a libraries <<< "${LAPACK_LIBS-$default_libraries}"
"${FC:-gfortran}" "${compile_flags[@]}" "${sources[@]}" \
   "$root/fortran_eos/app/nv_profile_generate.f90" "${link_flags[@]}" "${libraries[@]}" -o nv-profile
./nv-profile "$root/data/eos_original_papers_20261123/research_nv/nv1973_original_source.dat" "$tmp/tables" "${2:-128}" "${3:-49}"
cp -R tables/. "$out/"
