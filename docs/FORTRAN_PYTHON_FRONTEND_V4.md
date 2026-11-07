# Fortran EOS core and Python frontend v4

## Architecture

The default CompOSE data reader, EOSData serializer, CompOSE writer and C1-C11 thermodynamic verifier are implemented in Fortran under the fortran_eos directory. The Python package under python_frontend/compose_eos provides familiar public entrypoints and invokes the Fortran executable through a stable process, CSV and JSON contract. No Python thermodynamic fallback is silently used when the native compiler or executable is missing.

The TOV and relativistic Cowling interfaces live in python_frontend/compose_cowling/solver.py. They compile the existing Fortran case_driver and cell_solve sources in dependency order, then expose Python functions solve_tov_star, solve_tov_sequence and solve_cowling_spectrum. Python performs input validation, file preparation, process orchestration and result packaging; stellar structure and oscillation integration remain Fortran.

The native ABI is intentionally a file/CLI contract in this first port. This avoids compiler-specific derived-type ABI leakage and Fortran/Python array-order ambiguity. A C/ISO_C_BINDING or f2py layer can be added later without changing the user-facing Python objects.

## Current public API

    from compose_eos import EOSData, read_compose, write_compose, check_consistency, fill_density_gap
    from compose_cowling import solve_tov_star, solve_tov_sequence, solve_cowling_spectrum

    eos = read_compose("path/to/compose_table")
    report = check_consistency(eos)
    filled = fill_density_gap(eos, n_fill=60, gap_phase=7)  # pass gap_phase only for an explicitly approved phase crossing
    write_compose("path/to/output", eos)
    star = solve_tov_star(eos, pc_pa=1.5e34)
    sequence = solve_tov_sequence(eos, [1.0e34, 1.2e34, 1.5e34])
    modes = solve_cowling_spectrum(
        eos, pc_pa=1.4957e34,
        frequency_lo_hz=[2000, 5800], frequency_hi_hz=[3200, 8600]
    )

The first call builds the standalone EOS executable with gfortran if no explicit binary is supplied. The TOV/Cowling frontend similarly builds the existing Fortran solver once and caches it by source hash. COMPOSE_EOS_FORTRAN_EXE and COMPOSE_EOS_FORTRAN_ROOT can be set to pin the executable or source root.

## Thermodynamic representation

The CompOSE EOSData model retains original source n_B, Q1-Q7, optional Q-extra columns, phase labels and available composition fields. Its thermodynamic properties follow the original CompOSE definitions:
[
P=n_B Q_1,quad F=(Q_6+1)m_n,quad
mu_B=(Q_3+1)m_n,quad epsilon=n_B(Q_7+1)m_n.
]
The Fortran verifier checks source identities before exporting the normalized pressure-energy table. The stellar provider then consumes the v3 direct representation:
[
g(P)=epsilon(P)-P,quad u=ln(P/P_*),quad
c_s^2=rac{P}{P+g_u},quad
rac{dln n_B}{du}=rac{P+g_u}{2P+g}.
]
Original pressure and energy-density knots are retained; density and chemical potential are reconstructed by the first law. This interpolation closure is not a new microscopic EOS.

## Build and test

    python -m pip install -e .
    python -m pytest -q python_frontend/tests
    bash fortran_eos/build.sh build/compose_eos_fortran
    build/compose_eos_fortran read path/to/compose path/to/output
    build/compose_eos_fortran verify path/to/compose report.json

The standalone native build requires a Fortran 2018 compiler. The solver build currently targets gfortran and the included oscillation_source workspace. No FPM 0.13.0 dependency is required by the Python frontend.

## Important feature-parity boundary

The native CompOSE reader/writer and C1-C11 verifier are ported. The integral-preserving master-function density-gap fill is now implemented in Fortran and exposed as `fill_density_gap`; it retains source knots outside the gap and rejects an implicit cross-phase fill unless `gap_phase` is supplied. Full master-function smoothing/candidate selection, all six crust-core matching constructions, literature-specific EOS solvers, genuine Maxwell/Gibbs phase transitions, causal completion, and full TOV sequence diagnostics are not yet all ported to the new Fortran EOS core. Those remaining APIs raise `NativeFeatureUnavailable` rather than silently dispatching to the old Python implementation. The original Python package remains under eos_generator/compose_eos_pkg as a differential reference, not the installed backend.

The existing Fortran TOV/Cowling kernels are callable from Python, but solve_tov_sequence is currently a Python-controlled scan over independent Fortran TOV solves, and Cowling full-cohort convergence is still a separately audited campaign. This is a functional frontend, not a claim of complete one-to-one parity with every historic Python utility.
