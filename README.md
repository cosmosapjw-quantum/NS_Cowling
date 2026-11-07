# Integrated CompOSE EOS and relativistic Cowling pipeline v4

The default numerical path is Fortran. The Python frontend exposes CompOSE reading/writing, EOS thermodynamic verification, normalized EOS export, TOV solves and relativistic Cowling f/p1 calculations. See docs/FORTRAN_PYTHON_FRONTEND_V4.md and docs/FORTRAN_EOS_FEATURE_PARITY.md.

Install and test:

    python -m pip install -e .
    python -m pytest -q python_frontend/tests

Native EOS executable:

    bash fortran_eos/build.sh build/compose_eos_fortran
    build/compose_eos_fortran read path/to/compose path/to/internal
    build/compose_eos_fortran verify path/to/compose report.json

Python API:

    from compose_eos import read_compose, write_compose, check_consistency
    from compose_cowling import solve_tov_star, solve_tov_sequence, solve_cowling_spectrum

The legacy Python implementation remains in eos_generator/compose_eos_pkg/ only as a reference for differential tests. Native integral-preserving density-gap filling and six strict-thermo-only smooth crust-core bridges are implemented on the feature branch. Composition-aware crust/core reconstruction, literature_beta_eq, full master-function smoothing/candidate selection, literature-specific EOS solvers, genuine Maxwell/Gibbs phase-transition generation, and complete TOV sequence diagnostics remain out of scope. The legacy `maxwell` matching token is a smooth near-pressure plateau with no density jump. See docs/FORTRAN_EOS_CRUST_CORE_MATCHING_V1.md for equations, interface and limitations.

The stellar solver is nonrotating, spherical, perfect-fluid ell=2 polar modes in the relativistic Cowling approximation, with cold barotropic perturbations and a finite positive pressure surface cutoff. It is not a full-GR metric-perturbation calculation and does not compute gravitational-wave damping.


## Native smooth crust-core matching

The Fortran EOS frontend exposes the six historical smooth barotropic matching methods through `match_crust_core` in strict thermo-only mode. Inputs that contain composition require an explicit `drop_composition=True`; the native routine verifies cold beta-equilibrium, zero entropy, Q6=Q7, monotonicity, causality and the First-Law integral. The method token `maxwell` is a compatibility label for a smooth near-pressure plateau, not a true Maxwell density jump. See [the self-contained implementation note](docs/FORTRAN_EOS_CRUST_CORE_MATCHING_V1.md) and [the adversarial audit](docs/ADVERSARIAL_AUDIT_CRUST_CORE_MATCHING_20261107.md).
