#!/usr/bin/env bash
# Emit one CSV row per Fortran build, integration-control setting, star and mode.
set -euo pipefail
repo=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
inputs=${1:-"$repo/data/cowling_fortran_source_exact_root_c12_20261121"}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 LC_ALL=C

for build in checked optimized; do
    bash "$repo/scripts/build_c16_fortran_original_f_p1.sh" "$work/$build/c16" "$build" >&2
done

echo 'model,controls,build,mode,frequency_Hz,difference_from_driver_reference_Hz,source_knots'
for model in M2 M3; do
    case "$model" in
        M2)
            fref=1842.3116593617328
            half_step=1.06166266658578863e-3
            expected_sha=21fc726337dba746ba5d461c43a1c429138abd15fb66acb6349d2aeca4b4648e
            ;;
        M3)
            fref=1564.2448940763072
            half_step=7.39719202404320855e-4
            expected_sha=f85dcebf71dd36959a466328ef05e4fbf680cb3485cccd7eb06ddaf382f85dc4
            ;;
    esac
    source_file="$inputs/ORIGINAL_FORTRAN_ROOT_${model}.dat.gz"
    printf '%s  %s\n' "$expected_sha" "$source_file" | sha256sum --check >&2
    gzip -dc "$source_file" > "$work/original.dat"
    # Preserve every source and background record. The control record follows
    # the two EOS tables and background sizes in the C12 input format.
    awk -v step="$half_step" '
        NR == 1 { controls = $1 + $2 + 5 }
        NR == controls { $3 = "2e-10"; $4 = "2e-13"; $5 = step }
        { print }
    ' "$work/original.dat" > "$work/tighter.dat"
    for controls in original tighter; do
        for build in checked optimized; do
            "$work/$build/c16" "$work/$controls.dat" "$fref" 3 > "$work/$build.log"
            awk -v star="$model" -v controls="$controls" -v build="$build" '
                BEGIN { OFS = "," }
                $1 == "C16_FORTRAN_F_P1" { print star, controls, build, ($3 == 1 ? "f" : "p1"), $5, $7, $9 }
            ' "$work/$build.log"
        done
        cmp "$work/checked.log" "$work/optimized.log"
    done
done
