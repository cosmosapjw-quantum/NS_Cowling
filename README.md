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

The legacy Python implementation in eos_generator/compose_eos_pkg/ is preserved as a read-only characterization and differential reference; it is not a production calculation path. The target architecture is one Fortran-owned calculation backend with a Python frontend limited to input validation, data transport, native invocation, report decoding and presentation. Native smooth crust-core matching currently exposes css, cubic_F, cubic_P, quintic_F and quintic_P in thermo-only mode. The historical maxwell plateau alias is fail-closed because it has zero density jump. The separate sharp Maxwell kernel is currently restricted to composition-free branch inputs; full CompOSE composition, literature-derived beta-equilibrium, true Maxwell pair matching, and native TOV/tidal/Cowling integration into the EOS390 pipeline remain active migration tasks. See docs/EOS390_NATIVE_FORTRAN_UNIFICATION_CONTRACT_20261111.md.

The stellar solver is nonrotating, spherical, perfect-fluid ell=2 polar modes in the relativistic Cowling approximation, with cold barotropic perturbations and a finite positive pressure surface cutoff. It is not a full-GR metric-perturbation calculation and does not compute gravitational-wave damping.


## Native smooth crust-core matching

The native smooth matching API currently exposes five methods: css, cubic_F, cubic_P, quintic_F and quintic_P. The old Python `maxwell` plateau remains only in the read-only legacy package. The public native matcher rejects `method='maxwell'` until a branch-selected sharp transition with nonzero density jump and full CompOSE composition is routed into the same API. The native CompOSE writer now supports phase-only sidecars (`Npairs=0`, `Nquad=0`) so Iphase labels are not erased when abundance values are absent. See [the EOS390 native Fortran contract](docs/EOS390_NATIVE_FORTRAN_UNIFICATION_CONTRACT_20261111.md).


## Native first-order Maxwell transition

The separate native `maxwell-transition` kernel is a genuine one-conserved-baryon-charge construction for two locally neutral cold branches. Its current v1 scope rejects composition-bearing inputs and nonzero Nadd, so it is not yet the EOS390 production path. The repository contains phase-segmented Fortran TOV/Cowling solvers for sharp jumps, but their integration with CompOSE-generated Maxwell outputs must be verified on source EOS pairs. The EOS390 contract makes branch selection, composition-preserving transition output, and interface validation explicit prerequisites.


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
