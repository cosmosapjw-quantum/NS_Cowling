#!/usr/bin/env bash
# Independently generated 128/256 profile directories plus the built verifier.
set -euo pipefail
[[ $# == 3 ]] || { echo 'usage: script COARSE FINE NV_PROFILE_VERIFY' >&2; exit 2; }
root=$(cd "$(dirname "$0")/../.." && pwd)
coarse=$(realpath "$1")
fine=$(realpath "$2")
exe=$(realpath "$3")
testdir=$(mktemp -d)
trap 'rm -rf "$testdir"' EXIT
verify() { python3 "$root/scripts/verify_nv_reconstructed_profiles.py" "$1" "$2" --exe "$exe"; }
verify "$coarse" "$fine" >"$testdir/normal.json"
echo "PASS generated coarse/fine pair"
cp -R "$coarse" "$testdir/problem"
expect_failure() {
    if verify "$1" "$2" >"$testdir/output" 2>"$testdir/error"; then
        echo "Unexpected acceptance of $3" >&2
        exit 1
    fi
    grep -F -- "$3" "$testdir/error" >/dev/null || { cat "$testdir/error" >&2; exit 1; }
    echo "PASS rejected: $3"
}
corrupt_field() {
    local file=$1 field=$2 value=$3
    awk -F, -v name="$field" -v value="$value" 'BEGIN {OFS=","}
        NR==1 {for(i=1;i<=NF;i++) if($i==name) col=i; if(col==0) exit 2}
        NR==2 {$col=value} {print}' "$coarse/$file" >"$testdir/problem/$file"
}
corrupt_field uniform_beta.csv Ye 0.49
expect_failure "$testdir/problem" "$fine" 'uniform neutrality or nucleon closure'
corrupt_field uniform_beta.csv mu_p_MeV 800
expect_failure "$testdir/problem" "$fine" 'CSV beta affinity inconsistent'
corrupt_field uniform_beta.csv cs2_over_c2_FD 1.2
expect_failure "$testdir/problem" "$fine" 'uniform local sound speed'
corrupt_field uniform_beta.csv minimum_sampled_charge_curvature_MeV_fm3 -1
expect_failure "$testdir/problem" "$fine" 'uniform charge curvature'
cp "$coarse/uniform_beta.csv" "$testdir/problem/uniform_beta.csv"
corrupt_field printed.csv max_nucleon_density_fm3 NaN
expect_failure "$testdir/problem" "$fine" 'nonfinite numeric CSV field'
expect_failure "$coarse" "$coarse" 'Use independent coarse and fine profile runs'
if [[ -f "$coarse/profile_run.nml" && -f "$fine/profile_run.nml" ]]; then
    expect_failure "$fine" "$coarse" 'refined radial mesh must be finer'
fi
cp -R "$fine" "$testdir/fine_problem"
awk -F, 'BEGIN {OFS=","}
    NR==1 {for(i=1;i<=NF;i++) if($i=="E_minus_mn_MeV") col=i; if(col==0) exit 2}
    NR==2 {$col+=0.001} {print}' "$fine/printed.csv" >"$testdir/fine_problem/printed.csv"
expect_failure "$coarse" "$testdir/fine_problem" 'refined CompOSE energy per baryon'
