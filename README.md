# Integrated CompOSE EOS and relativistic Cowling pipeline v4

The default numerical path is Fortran. The Python frontend exposes CompOSE reading/writing, EOS thermodynamic verification, normalized EOS export, TOV solves and relativistic Cowling f/p1 calculations. See docs/FORTRAN_PYTHON_FRONTEND_V4.md and docs/FORTRAN_EOS_FEATURE_PARITY.md.

The 2026-11-23 EOS generator and matching re-audit, expanded all-source campaign, numerical scope, and delivery instructions are documented in [the current audit report](docs/EOS_REAUDIT_ALL_PAIRS_20261123.md). Earlier native generation work is documented in [the 2026-11-17 research report](docs/research_matching_20261117/REPORT_KO.md), [the 18-model closure matrix](data/research_matching_20261117/generator/RETAINED_COHORT_MATRIX.json), and [reproduction and handoff instructions](docs/research_matching_20261117/HANDOFF_KO.md). Production model evaluation, minimization and matching run in Fortran; Python validates requests, transports data and decodes reports. Independent research oracles and the legacy Python implementation remain explicitly separate. Original-source fidelity, independent physical validation, numerical support and CompOSE comparison are separate judgments; exact CompOSE agreement is not an unconditional pass criterion.

Install and test:

    python -m pip install -e .
    python -m pytest -q python_frontend/tests

Native EOS executable:

    bash fortran_eos/build.sh build/compose_eos_fortran
    build/compose_eos_fortran/compose_eos_fortran read path/to/compose path/to/internal
    build/compose_eos_fortran/compose_eos_fortran verify path/to/compose report.json

Source-checkout EOS generation (fresh output directories):

    PYTHONPATH=python_frontend python3 -m compose_eos.cli generators
    PYTHONPATH=python_frontend python3 -m compose_eos.cli generate GM_GM1 build/gm1_run --points 161
    PYTHONPATH=python_frontend python3 -m compose_eos.cli generate SPG_M1 build/spg_m1_run --points 257 --variant AME2408_candidate
    PYTHONPATH=python_frontend python3 -m compose_eos.cli generate SPG_M5 build/spg_m5_inferred --points 129 --variant public_composition_inferred_20261116

The catalog exposes five core models, five SPG candidate inner-crust models and eight GPPVA spatial-TF models. TF generation requires explicit density and radial-resolution grids. Generation success is limited to requested rows and their native output verification; it does not certify original SPG calibration, a global minimum, a complete joined crust, or arbitrary interpolation across phase boundaries. Density-dependent physical support and composition availability accompany the outputs. The conditional inferred SPG variant uses a versioned public-composition reconstruction; the default remains AME2408_candidate, and neither variant asserts the missing original calibration. Explicit unsupported phase/segment boundaries are enforced when matching. All 18 declared native generation routes have execution evidence; the overall original-model fidelity and continuous full-crust goal remains **PARTIALLY_COMPLETE**.

Python API:

    from compose_eos import read_compose, write_compose, check_consistency
    from compose_cowling import solve_tov_star, solve_tov_sequence, solve_cowling_spectrum

The legacy Python implementation in eos_generator/compose_eos_pkg/ is preserved as a read-only characterization and differential reference; it is not a production calculation path. Production EOS generation, thermodynamic conversion, matching and physical verification run in Fortran. Python handles requests, file transport, native invocation, report aggregation and presentation. The seven canonical matching methods and their exact contracts are described in [the current matching specification](docs/FORTRAN_EOS_CRUST_CORE_MATCHING_V1.md). The native API includes a genuine sharp Maxwell transition; it does not route that name to the historical zero-jump pressure plateau. Matching views preserve physical columns and phase labels, and any removal of auxiliary columns or abundances is explicit and recorded. No abundances are fabricated.

The stellar solver is nonrotating, spherical, perfect-fluid ell=2 polar modes in the relativistic Cowling approximation, with cold barotropic perturbations and a finite positive pressure surface cutoff. It is not a full-GR metric-perturbation calculation and does not compute gravitational-wave damping.


## Native crust-core matching

The native API exposes `cubic_F`, `quintic_F`, `cubic_P`, `quintic_P`, `first_law_constrained_hermite_pressure_bridge`, `css`, and `maxwell`. The first five construct continuous bridges; CSS can carry a crust-side density jump, and Maxwell solves equality of pressure and chemical potential with a nonzero density jump. The old Python pressure-plateau alias remains confined to the legacy package.

The 2026-11-23 matcher retries rejected construction windows and exhausts the supplied eligible endpoint grid before declaring that search unsuccessful. Polynomial pieces are checked at their floating-point extrema; the first-law pressure bridge uses its analytic interval certificate as well as diagnostic samples. These are floating-point certificates, not rigorous interval-arithmetic proofs. Retained source domains and serialized outputs receive separate checks. Failed construction, unsupported multiple coexistence roots, missing source inputs, and post-construction admission failures remain distinct outcomes.

## Native first-order Maxwell transition

The native Maxwell kernel constructs coexistence between two locally neutral cold branches with one conserved baryon charge. Its thermodynamic input contract excludes abundance columns and nonzero Nadd; copied matching views may explicitly project these columns while retaining the original source files and phase labels. The CompOSE writer preserves phase-only sidecars (`Npairs=0`, `Nquad=0`). Each admitted campaign output is checked again after writing and reading CompOSE, including independent interface admission. A selected source endpoint chain is not a certification of a zero-pressure surface, every possible interpolation, a global nuclear ground state, or a stable stellar sequence. Finite phase-component matches are listed separately from selected-source endpoint chains. See [the physical audit](docs/EOS_PHYSICS_SOURCE_AUDIT_20261123.md) and [the current campaign report](docs/EOS_REAUDIT_ALL_PAIRS_20261123.md).

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

## Sharp Maxwell phase-interface TOV/Cowling reference path (2026-11-09)

The sharp-interface follow-up is implemented on the `research/maxwell-reference-validation-v1` branch and tracked in [PR #2](https://github.com/cosmosapjw-quantum/NS_Cowling/pull/2). It consumes the existing native one-charge Maxwell construction as two separate phase providers, integrates the high-density core to exact `P_t`, and restarts the low-density envelope without interpolating across the density jump. Its `ell=2` Cowling reference path distinguishes slow conversion (`[xi^r]=[Delta P]=0`) from rapid conversion (`[Delta P]=0`, `[xi^r-Delta P/P'_0]=0`).

**Implementation boundary (updated 2026-11-09):** the phase-segmented stellar TOV/Cowling reference solver remains available in Python/SciPy, and native Fortran now exposes production generic `solve_tov_star` and `solve_cowling_spectrum` entry points for both smooth and one-sharp-Maxwell backgrounds. The native Maxwell route accepts two already-coexisting phase EOS providers and preserves separate core/envelope backgrounds. Its Cowling spectrum route performs a bounded frequency-grid scan and bracket refinement for slow or rapid conversion; it does not yet assign `f`, `p1`, or interface-mode identities. No realistic CompOSE hybrid-EOS result, stable-branch, maximum-mass, twin-star or full-GR claim is made. The strong-jump analytic fixture is for interface-condition testing only; its rapid-onset Seidov diagnostic is unstable. The small-jump fixture validates convergence toward the smooth-EOS solution.

See [the interface derivation](docs/MAXWELL_TOV_COWLING_INTERFACE_V1.md), [the independent validation record](docs/MAXWELL_REFERENCE_VALIDATION_V1_20261109.md), and [the final adversarial audit](docs/ADVERSARIAL_AUDIT_MAXWELL_TOV_COWLING_20261109.md).

## Native physical EOS generation (2026-11-16)

The public generator now runs core5, SPG5 and spatial TF8 physics in Fortran, with a Python frontend. See [native regeneration and limits](docs/research_native_20261116/REPRODUCE.md), [scientific report](docs/research_native_20261116/REPORT_KO.md), and [local handoff](docs/research_native_20261116/HANDOFF_KO.md). Native execution for all18 models does not certify unresolved author surface fits, full microscopic BHF/ETFSI reproduction, or continuous mesh-converged crust support.
