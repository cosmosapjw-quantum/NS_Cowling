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

**Implementation boundary (updated 2026-11-09):** the phase-segmented stellar TOV/Cowling reference solver remains available in Python/SciPy, and native Fortran now exposes production generic `solve_tov_star` and `solve_cowling_spectrum` entry points for both smooth and one-sharp-Maxwell backgrounds. The native Maxwell route accepts two already-coexisting phase EOS providers and preserves separate core/envelope backgrounds. Its Cowling spectrum route performs a bounded frequency-grid scan and bracket refinement for slow or rapid conversion; it does not yet assign `f`, `p1`, or interface-mode identities. No realistic CompOSE hybrid-EOS result, stable-branch, maximum-mass, twin-star or full-GR claim is made. The strong-jump analytic fixture is for interface-condition testing only; its rapid-onset Seidov diagnostic is unstable. The small-jump fixture validates convergence toward the smooth-EOS solution.

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



The native sharp-Maxwell Fortran path is tracked independently in [the implementation contract](docs/MAXWELL_NATIVE_FORTRAN_V1.md), [the adversarial audit](docs/ADVERSARIAL_AUDIT_MAXWELL_NATIVE_FORTRAN_20261109.md), and [the research DAG](docs/MAXWELL_NATIVE_FORTRAN_DAG_20261109.md). Typed production generics now dispatch smooth and one-interface Maxwell TOV/Cowling paths. The latest exact-head CI checks 12 cell, 22 analytic Maxwell, and 74 normalized-table assertions; builds the production driver; tests slow and rapid spectra with two positive roots each; varies ODE and root/grid settings; and compares native TOV/Cowling roots against the independent Python/SciPy phase-aware reference. The generated analytic fixture only is validated: root-family identity remains unresolved and no source-backed realistic CompOSE hybrid EOS has been run. PR #3 remains draft/open/unmerged.

## Native production solver dispatch (2026-11-09)

The Fortran module `ns_production_solvers` exposes generic production entry points:
- `solve_tov_star(eos, pc, options, profile [, eos_events])` preserves the existing smooth-EOS integration and optional smooth pressure-stop contract.
- `solve_tov_star(eos_high, eos_low, pc, pt, options, maxwell_star)` invokes the phase-segmented sharp-Maxwell TOV construction.
- `solve_cowling_spectrum(background, grid, options, spectrum)` and the corresponding `star_profile_t` overload preserve the smooth path.
- `solve_cowling_spectrum(maxwell_star, grid, options, conversion, spectrum)` scans a two-phase background separately for `slow` or `rapid` conversion.

The Maxwell scanner retains trial statuses, refines sign/near brackets with bisection, reruns Cowling shots at roots, and checks final surface residuals. Exact-head CI now verifies two slow and two rapid roots on the generated analytic fixture. Separate ODE refinement moves roots by at most 3.2810e-5 relative; denser-grid/tighter-bisection refinement moves them by at most 2.6561e-10. The independent Python/SciPy reference agrees in TOV M/R/interface radius at 5.5e-10 or better and in the four unlabelled Cowling roots within 5.21e-5 relative. The roots remain family-unresolved; the finite grid is not a spectrum-completeness proof. See [the DAG](docs/MAXWELL_NATIVE_FORTRAN_DAG_20261109.md), [the adversarial audit](docs/ADVERSARIAL_AUDIT_MAXWELL_NATIVE_FORTRAN_20261109.md), and [the local-agent resume prompt](docs/LOCAL_AGENT_RESUME_MAXWELL_PRODUCTION_DISPATCH_20261109.md).


The latest exact-head workflow evidence is recorded in the adversarial audit and local-agent resume file. The successful push and pull-request runs both check the workflow/docs head while all native Fortran source code remains at the separately verified implementation commit `f8678cb13a18f63bad9f8db65c4cd3e66e491ad3`.


## Native sharp-Maxwell production API and positive-root gate (2026-11-09)

The Fortran production generic APIs now dispatch both smooth and one-interface sharp-Maxwell TOV/Cowling paths. Exact-head CI uses generated analytic phase tables and verifies two bisection-refined roots for each conversion prescription. The numerical roots are:
- Slow: 1764.8989380477 Hz and 4937.6449322607 Hz.
- Rapid: 4937.0162338857 Hz and 19493.6326018255 Hz.

These values are synthetic regression-fixture roots only. All four mode-family labels remain `unresolved`; root convergence, Python/SciPy cross-validation, realistic CompOSE phases and mode identity remain open. The DAG and adversarial audit contain exact CI provenance.


## Native Maxwell convergence and independent reference gate (2026-11-09)

Exact tested head: `afcdf948f0ca6844a6c3402f66e184a46aa78978` (Fortran/table-test increment `9d3ef4bfe8b018b17887504d60e25d1c94971b6e`). Both push and PR Actions passed. The generated two-phase fixture now has explicit ODE-refinement and frequency-grid/root-tolerance-refinement runs for slow and rapid conversion, plus an independent Python/SciPy phase-aware cross-check.

- Maximum root shift under the tighter ODE setting: (3.2810\times10^{-5}).
- Maximum root shift under the denser grid and tighter bisection setting: (2.6561\times10^{-10}).
- Python/SciPy TOV relative differences: mass (5.47\times10^{-10}), radius (2.15\times10^{-10}), interface radius (2.67\times10^{-10}).
- Maximum Fortran/Python Cowling root-set relative difference: (5.205\times10^{-5}).

All four roots remain `unresolved` by mode family. This is a numerical verification on generated analytic tables, not a realistic neutron-star spectrum or source-backed CompOSE result. Exact artifacts and acceptance limits are recorded in the DAG and adversarial audit.


## Mass-dependent smooth crust-core matching spectra (2026-11-09)

A separate follow-up track is maintained on research/mass-dependent-cowling-v1. It extends the frozen 234-row smooth-matching cohort across the stellar mass sequence rather than reusing the old single-mass spectrum. The 234 rows are 39 fixed core/crust pairs crossed with six matching methods, not 234 independent EOS samples. The historical method token maxwell remains a smooth pressure-plateau bridge, not a sharp first-order density jump.

The new workflow scans TOV central pressure, extracts the initial connected rising-mass branch, rejects gaps and out-of-branch target masses, then sends adjacent stable-branch pressure brackets to the native mass-root solver before calculating fresh Cowling spectra. Six-method paired deviations are admitted only for complete core/crust/mass cells. The predeclared mass grid is 0.8, 1.0, 1.2, 1.33, 1.4, 1.6, 1.8 and 2.0 solar masses.

See the research plan at docs/MASS_DEPENDENT_COWLING_RESEARCH_PLAN_20261109.md, DAG at docs/MASS_DEPENDENT_COWLING_DAG_20261109.md, and local-agent resume prompt at docs/LOCAL_AGENT_RESUME_MASS_DEPENDENT_COWLING_20261109.md. The pilot/full GitHub Actions workflow is .github/workflows/mass-dependent-cowling.yml. This is an implementation increment awaiting exact-head CI and artifact review; no new mass-dependent scientific result or maximum-mass claim is made yet.
