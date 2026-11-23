#!/usr/bin/env bash
# Replay pre-frozen catalogs calling two distinct native executable builds.
set -euo pipefail
if (( $# < 5 )); then
    echo "Usage: $0 BASELINE_BINARY CANDIDATE_BINARY CATALOG_DIR NEW_SCRATCH_DIR LANE..." >&2
    exit 2
fi
root=$(cd "$(dirname "$0")/../.." && pwd)
baseline=$(realpath "$1")
candidate=$(realpath "$2")
catalogs=$(realpath "$3")
runtime=$(realpath -m "$4")
shift 4
if [[ -e "$runtime" || "$runtime" == "$root" || "$runtime" == "$root/"* ]]; then
    echo "Supply a fresh runtime directory outside the repository." >&2
    exit 2
fi
umask 077
mkdir -p "$runtime"
for arm in baseline candidate; do
    exe=$baseline
    if [[ $arm == candidate ]]; then exe=$candidate; fi
    for lane in "$@"; do
        python3 "$root/scripts/run_eos_all_pairs_20261123.py" \
            --catalog "$catalogs/$lane.json" --native-exe "$exe" \
            --out "$runtime/$arm/run/$lane" --private-artifacts "$runtime/$arm/private/$lane" \
            --work "$runtime/$arm/tmp" --n-bridge 80 --workers "${REPLAY_WORKERS:-3}"
    done
done
python3 "$root/fortran_eos/tests/compare_matching_replays.py" \
    "$runtime/baseline/run" "$runtime/candidate/run" "$runtime/comparison.json" "$@"
