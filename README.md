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

The legacy Python implementation remains in eos_generator/compose_eos_pkg/ only as a reference for differential tests. Native integral-preserving density-gap filling and six strict-thermo-only smooth crust-core bridges are implemented on the feature branch. A native cold master-function reconstruction v1 is now exposed through the Python frontend. It uses local weighted-cubic/PCHIP regularisation, reconstructs electron/muon equilibrium only when composition is present, and applies a bounded electron-plus-lattice correction on eligible supplied outer-crust rows. This is a partial port, not the historical global-spline/candidate-selection policy or a microscopic nuclear solver. A separate native one-conserved-charge Maxwell construction is available for two caller-confirmed locally neutral cold phases, but remains thermo-only and is not verified for phase-interface TOV sequences. Composition-aware crust-core matching, literature_beta_eq, Gibbs mixed-phase generation, Seidov onset evaluation and phase-interface/TOV stability remain outside the verified native path. The legacy `maxwell` matching token remains a smooth near-pressure plateau with no density jump. See docs/FORTRAN_EOS_CRUST_CORE_MATCHING_V1.md and docs/FORTRAN_EOS_MAXWELL_TRANSITION_V1.md.

The stellar solver is nonrotating, spherical, perfect-fluid ell=2 polar modes in the relativistic Cowling approximation, with cold barotropic perturbations and a finite positive pressure surface cutoff. It is not a full-GR metric-perturbation calculation and does not compute gravitational-wave damping.


## Native smooth crust-core matching

The Fortran EOS frontend exposes the six historical smooth barotropic matching methods through `match_crust_core` in strict thermo-only mode. Inputs that contain composition require an explicit `drop_composition=True`; the native routine verifies cold beta-equilibrium, zero entropy, Q6=Q7, monotonicity, causality and the First-Law integral. The method token `maxwell` is a compatibility label for a smooth near-pressure plateau, not a true Maxwell density jump. See [the self-contained implementation note](docs/FORTRAN_EOS_CRUST_CORE_MATCHING_V1.md) and [the adversarial audit](docs/ADVERSARIAL_AUDIT_CRUST_CORE_MATCHING_20261107.md).


## Native first-order Maxwell transition

The feature branch adds `construct_maxwell_transition` and the native `maxwell-transition` CLI for a single globally conserved baryon charge with local electric neutrality imposed separately in each supplied phase. The caller must explicitly confirm local charge neutrality because a composition-free barotrope cannot prove it. The Fortran kernel uses a First-Law-consistent cubic-Hermite representation in baryon chemical potential, analytic whole-interval mechanical and causal checks, and stationary-point-isolated pressure-root searches. It rejects composition-bearing inputs, nonzero Nadd, and mismatched charge-fraction slices. Its output is an internal phase-labelled barotrope, not a composition-preserving CompOSE product; both Python and native Fortran writers refuse a complete CompOSE export without phase composition. Gibbs, Seidov stability evaluation and phase-interface TOV support remain unported or unverified. See [the Maxwell transition note](docs/FORTRAN_EOS_MAXWELL_TRANSITION_V1.md) and [the adversarial audit](docs/ADVERSARIAL_AUDIT_MAXWELL_TRANSITION_20261107.md).


## Sharp Maxwell phase-interface TOV/Cowling reference path (2026-11-09)

The sharp-interface follow-up is implemented on the `research/maxwell-reference-validation-v1` branch and tracked in [PR #2](https://github.com/cosmosapjw-quantum/NS_Cowling/pull/2). It consumes the existing native one-charge Maxwell construction as two separate phase providers, integrates the high-density core to exact `P_t`, and restarts the low-density envelope without interpolating across the density jump. Its `ell=2` Cowling reference path distinguishes slow conversion (`[xi^r]=[Delta P]=0`) from rapid conversion (`[Delta P]=0`, `[xi^r-Delta P/P'_0]=0`).

**Implementation boundary:** the phase-segmented stellar TOV/Cowling reference solver is Python/SciPy. Native Fortran now independently implements and tests both interface-transfer maps, but the production Fortran `solve_tov_star` and `solve_cowling_spectrum` remain smooth-EOS-only. No stable-branch, maximum-mass, twin-star or full-GR claim is made. The strong-jump analytic fixture is for interface-condition testing only; its rapid-onset Seidov diagnostic is unstable. The small-jump fixture validates convergence toward the smooth-EOS solution.

See [the interface derivation](docs/MAXWELL_TOV_COWLING_INTERFACE_V1.md), [the independent validation record](docs/MAXWELL_REFERENCE_VALIDATION_V1_20261109.md), and [the final adversarial audit](docs/ADVERSARIAL_AUDIT_MAXWELL_TOV_COWLING_20261109.md).

## Native master-function reconstruction v1

    from compose_eos import reconstruct_master
    rebuilt, report = reconstruct_master(
        eos, smooth=True, smooth_sigma=0.05,
        fix_dissolution=True, dissolution_tol=0.015,
        recompute_extra_cs2=False, return_report=True,
    )

The Fortran path defines \(F=f/n_B\) on \(x=\ln n_B\) and derives \(Q_1=F'\), \(\mu_H=F+Q_1=d\epsilon/dn_B\), \(P=n_BQ_1\), and \(\epsilon=n_BF\) from the same interpolant. The outer-crust bridge is now checked over the complete represented polynomial interval using stationary-point isolation; the former fixed mesh has been removed. The native input contract additionally enforces cold Q6=Q7 and rowwise charge-aware CompOSE Q3/Euler closure, including the Il=0 charge-potential convention; the barotropic sound-speed denominator remains mu_H=d epsilon/dn_B. Composition-free inputs do not receive fabricated particle fractions. See [the self-contained derivation and limits](docs/FORTRAN_EOS_MASTER_RECONSTRUCTION_V1.md) and [the bounded adversarial audit](docs/ADVERSARIAL_AUDIT_MASTER_RECONSTRUCTION_20261108.md).


## v4 physical and verification guide

See [the self-contained Fortran physics guide](docs/NS_COWLING_PHYSICS_AND_CODE_GUIDE_V4.md), [the final adversarial audit](docs/ADVERSARIAL_AUDIT_V4_FINAL_20261108.md) and [the machine-readable trace](docs/ADVERSARIAL_AUDIT_TRACE_V4_20261108.json).


## Verified v4 Fortran snapshot

Verified code head: \`1dbb75c27b704209cfeef43f7c49be0cfb28fd44\`. GitHub Actions reports 67 Python tests passed, the native Fortran polynomial-extrema test passed, all 69 tracked SHA-256 entries passed, and the source archive was uploaded. See the [final adversarial audit](docs/ADVERSARIAL_AUDIT_V4_FINAL_20261108.md) for run IDs and artifact digests.


## Native Fortran sharp-Maxwell TOV/Cowling v1

A separate phase-segmented native Fortran path now accepts independent high-density and low-density EOS tables, enforces the Maxwell pressure/chemical-potential/first-law conditions at the sharp interface, integrates TOV on two non-crossing pressure segments, and performs a single-frequency Cowling shot with slow or rapid phase-conversion transfer. See [the implementation and validation contract](docs/MAXWELL_NATIVE_FORTRAN_V1.md). This path is separate from the legacy smooth maxwell bridge and the 234-row v5 smooth-matching cohort. Its analytic fixture and exact-head CI gate validate the bounded implementation scope only; they do not yet establish a realistic hybrid-EOS spectrum survey.
