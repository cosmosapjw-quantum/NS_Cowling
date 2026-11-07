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

The native CompOSE reader/writer and C1-C11 verifier are ported. The integral-preserving master-function density-gap fill is implemented in Fortran and exposed as `fill_density_gap`; it retains source knots outside the gap and rejects an implicit cross-phase fill unless `gap_phase` is supplied. The six historical smooth crust-core bridges are implemented in Fortran and exposed as `match_crust_core` in strict thermo-only mode. This wrapper retains the legacy `(EOSData, info)` return contract and includes the native bridge report and independent C1--C11 consistency record. If either input carries composition, callers must explicitly pass `drop_composition=True`. The historic `maxwell` token means a smooth near-pressure plateau with no density jump, not a physical Maxwell coexistence construction. A separate one-conserved-charge Maxwell transition is now available through `construct_maxwell_transition`; it requires cold beta-equilibrated single-phase inputs and an explicit caller assertion that both phases are locally charge neutral. It returns a phase-labelled barotrope but does not preserve composition, evaluate the Seidov criterion, or certify phase-interface TOV integration. Python and native Fortran CompOSE writers reject a complete export without phase composition. Gibbs mixed phases, composition-aware interpolation, the `literature_beta_eq` mode, full master-function smoothing/candidate selection and literature-specific EOS solvers remain unported. Those APIs remain unavailable rather than silently dispatching to the old Python implementation. The original Python package remains under eos_generator/compose_eos_pkg as a differential reference, not the installed backend.

The existing Fortran TOV/Cowling kernels are callable from Python, but solve_tov_sequence is currently a Python-controlled scan over independent Fortran TOV solves, and Cowling full-cohort convergence is still a separately audited campaign. This is a functional frontend, not a claim of complete one-to-one parity with every historic Python utility.


## Native crust-core matching v1

See [FORTRAN_EOS_CRUST_CORE_MATCHING_V1.md](FORTRAN_EOS_CRUST_CORE_MATCHING_V1.md) for the standalone thermodynamic derivation, transition-window selector, six bridge definitions, native CLI, Python API, validation gates, benchmark provenance and limitations. The feature-parity matrix is updated conservatively: the matching port remains PARTIAL until the matching-specific differential CI has passed on the exact branch head.


## Native one-charge Maxwell transition v1

    from compose_eos import construct_maxwell_transition
    hybrid, report = construct_maxwell_transition(
        hadronic_phase, high_density_phase,
        assume_local_charge_neutrality=True, return_report=True,
    )

The Fortran kernel interpolates (P(\mu_B)) with endpoint slopes (n_B=dP/d\mu_B), reconstructs (epsilon=\mu_B n_B-P), enforces mechanical stability and causality over whole Hermite intervals, isolates all pressure-equality roots through stationary points of the cubic pressure difference, and constructs two equal-(P)/equal-(\mu_B) rows with distinct densities. The caller assertion is necessary because a composition-free table cannot prove local neutrality. Inputs with composition, nonzero Nadd or inconsistent charge-fraction metadata are rejected. The output is an internal phase-labelled barotrope, not a complete CompOSE release. See [FORTRAN_EOS_MAXWELL_TRANSITION_V1.md](FORTRAN_EOS_MAXWELL_TRANSITION_V1.md) and [ADVERSARIAL_AUDIT_MAXWELL_TRANSITION_20261107.md](ADVERSARIAL_AUDIT_MAXWELL_TRANSITION_20261107.md) for derivation, tolerances and limitations.
