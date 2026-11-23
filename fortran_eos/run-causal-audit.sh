#!/usr/bin/env bash
# Read physical cold-EOS CSVs; write one summary per input to stdout.
set -eu
dir=$(cd "$(dirname "$0")" && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
"${FC:-gfortran}" -std=f2018 -O2 -ffree-line-length-none "$dir/tools/cold_causal_endpoint_audit.f90" -o "$tmp/cold_causal_endpoint_audit"
"$tmp/cold_causal_endpoint_audit" "$@"
