#!/usr/bin/env bash
# Isolated diagnostic-only copy of the adapter; production source unchanged.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD="${1:?isolated build directory required}"
mkdir -p "$BUILD"
python3 - "$ROOT" "$BUILD" <<'PYBUILD'
import sys
from pathlib import Path
root,build=map(Path,sys.argv[1:])
s=(root/'fortran_eos/src/compose_eos_optimize_native.f90').read_text()
needle='    call diagnostics()\n    select case(mode)'
assert s.count(needle)==1
s=s.replace(needle,"""    call diagnostics()
    print '(a,4(1x,i0),11(1x,es25.16e3))','SLSQP_DIAGNOSTICS',mode,iter,ireset,line, &
      result%max_constraint_error,result%max_bound_violation,f,f0,gs,h1,h2,h3,tol,alpha,acc
    select case(mode)""")
(build/'compose_eos_optimize_native_diagnostic.f90').write_text(s)
PYBUILD
gfortran -std=legacy -O2 -g -ffixed-line-length-none -c "$ROOT/fortran_eos/vendor/slsqp/slsqp_optmz.f" -o "$BUILD/slsqp.o"
for module in compose_eos_types compose_eos_optimize_native compose_eos_tf_models compose_eos_spatial_tf; do
 source="$ROOT/fortran_eos/src/$module.f90"
 if [[ "$module" == compose_eos_optimize_native ]]; then source="$BUILD/compose_eos_optimize_native_diagnostic.f90"; fi
 gfortran -std=f2018 -O2 -g -ffree-line-length-none -fcheck=all -J"$BUILD" -I"$BUILD" -c "$source" -o "$BUILD/$module.o"
done
gfortran -std=f2018 -O2 -g -ffree-line-length-none -fcheck=all -I"$BUILD" "$ROOT/fortran_eos/tests/tf_finish_single_state.f90" \
 "$BUILD/compose_eos_types.o" "$BUILD/compose_eos_optimize_native.o" "$BUILD/compose_eos_tf_models.o" "$BUILD/compose_eos_spatial_tf.o" "$BUILD/slsqp.o" -llapack -lblas -o "$BUILD/tf_finish_single_state"
