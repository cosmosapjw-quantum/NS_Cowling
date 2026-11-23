#!/usr/bin/env bash
# Run with a freshly built native profile executable as argument 1.
set -euo pipefail
exe=$(realpath "${1:?usage: nv_profile_parallel_equivalence.sh NV_PROFILE_EXECUTABLE}")
root=$(cd "$(dirname "$0")/../.." && pwd)
source_file="$root/data/eos_original_papers_20261123/research_nv/nv1973_original_source.dat"
out=$(mktemp -d)
trap 'rm -rf "$out"' EXIT

for threads in 1 4; do
    mkdir -p "$out/$threads"/{printed,d03x100,uniform_beta}
    OMP_DYNAMIC=FALSE OMP_NUM_THREADS=$threads "$exe" "$source_file" "$out/$threads" 64 11
done
diff -qr "$out/1" "$out/4"
