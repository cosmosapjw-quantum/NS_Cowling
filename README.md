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

The legacy Python implementation remains in eos_generator/compose_eos_pkg/ only as a reference for differential tests. Full generator feature parity is not claimed: master-function reconstruction, all six crust-core matching methods, literature-specific EOS models, true Maxwell/Gibbs phase-transition generation and causal completion remain explicitly marked NOT_PORTED or PARTIAL.

The stellar solver is nonrotating, spherical, perfect-fluid ell=2 polar modes in the relativistic Cowling approximation, with cold barotropic perturbations and a finite positive pressure surface cutoff. It is not a full-GR metric-perturbation calculation and does not compute gravitational-wave damping.
