# Installed wheel verification

A fresh venv installed the local wheel with pip. From a directory outside the source checkout, the Python frontend located the installed Fortran sources under the environment prefix and built its own binary in an initially empty user cache. The venv inherited the environment's NumPy and SciPy installations. The host provided GNU Fortran 14.2 and the versioned LAPACK and BLAS shared libraries. The build used the default checked -O2 flags and a single thread.

The official Python console script read 1,028 rows from public reconstructed BPS + NV + D1MS cell 000001 and listed 29 generator representations. Both layered solves through the installed Python API agreed exactly with all 13 numeric columns checked in the stored standalone campaign. The isolated calls supplied no independent endpoint audit, so the causal status in the receipt is `unaudited`. The public campaign records the first failing interval separately.

The installed-package check exercises the public source archive distributed with the repository. Run from a full checkout with setuptools, wheel, NumPy, SciPy, gfortran and the LAPACK/BLAS development libraries available.

```bash
set -eu
root=$PWD
tmp=$(mktemp -d)
python -m pip wheel --no-deps --no-build-isolation --no-index -w "$tmp/dist" "$root"
python -m venv --system-site-packages "$tmp/venv"
"$tmp/venv/bin/python" -m pip install --no-deps --no-index "$tmp"/dist/*.whl
cd "$tmp"
unset PYTHONPATH COMPOSE_EOS_FORTRAN_EXE COMPOSE_EOS_FORTRAN_ROOT COMPOSE_EOS_NATIVE_BUILD_DIR
export XDG_CACHE_HOME="$tmp/cache" OMP_NUM_THREADS=1
src="$root/data/eos_reconstructed_matching_20261124/remaining_matches/BPS_MTW1964_MS1966_NV1973_Hermite/PUBLIC_MATCHED_EOS.zip"
ref="$root/data/eos_native_upgrade_20261124/layered_stars/public_stars.csv"
"$tmp/venv/bin/python" -I "$root/fortran_eos/tests/smoke_installed_layered.py" "$tmp"/dist/*.whl "$src" "$ref"
```

On the verification host the linker override was `NATIVE_LIBS="-l:liblapack.so.3 -l:libblas.so.3"`. WHEEL_SMOKE.json records this build's wheel hash, source data hashes, flags, numerical outputs and comparison coverage. Wheel ZIP timestamps may differ across builds.
