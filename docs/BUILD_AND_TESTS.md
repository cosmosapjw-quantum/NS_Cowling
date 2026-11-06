# v3 build and test instructions

## Python

Install the bundled generator and project:

    python -m pip install -e eos_generator/compose_eos_pkg
    python -m pip install -e .

Run the full project pipeline tests:

    python -m pytest -q pipeline/tests

Audit the original EOS tables and write non-mutating receipts:

    python scripts/compose_cowling.py audit-direct --source-dir data/eos_original_normalized --output-dir data/eos_direct_provider_v3

## Fortran

The default workspace is `oscillation_source/solutions`, with `oscillation_source/support` as its dependency. The repository's FPM build contract pins FPM 0.13.0. This target machine did not have FPM installed, so the official FPM harness was unavailable. The v3 source was instead compiled independently with gfortran using Fortran 2018, bounds/runtime checks, module dependency order, and a full link of both `case_driver` and `cell_solve`. This manual build is recorded separately and must not be reported as an FPM test pass.

Fortran direct-provider cohort audit (after building `case_driver`):

    python scripts/audit_fortran_direct_provider.py --driver /path/to/case_driver --source-dir data/eos_original_normalized --output-dir data/eos_direct_provider_v3 --bad-fixture oscillation_source/data/normalized/compose_134_si.csv

## Existing spectrum

The old q-linear mode cohort can be checked without recomputing it:

    python scripts/compose_cowling.py verify-existing --release oscillation_release_v1

This command does not make that legacy result set a v3 spectrum.
