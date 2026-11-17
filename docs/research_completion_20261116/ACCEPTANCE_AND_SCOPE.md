# Independent EOS regeneration: acceptance and model scope

Base c8ff2d6c97a93f1ccf621051861d3091bfe94d81; user instruction 2026-11-16.
The original objective is usable independent generation of five core and thirteen crust models, with composition. Core–crust matching and stellar observables are outside this loop.

## Three independent questions

1. **Physical implementation:** Does the executable evaluate the declared published energy functional, constraints, species, boundary conditions and phase prescription? Check analytic/numerical derivatives, conservation, thermodynamic identities, limiting cases and independent reference calculations.
2. **Numerical support:** On which density, resolution and model-parameter domain do actual calculations converge? Preserve optimizer failures, boundary/KKT states, mesh and density refinement, alternate seeds and phase competitors. A finite supported interval is usable within its recorded assumptions; an acausal high-density extension does not reject its lower-density part.
3. **Reference comparison:** What differs from CompOSE or another supplied table? Record identity, thermodynamic/composition conventions and deviations. Exact numerical agreement is not an unconditional admission criterion. A discrepancy is a diagnostic that may reflect an implementation defect, a different model/phase convention, rounded inputs, interpolation or a defect in the reference itself.

Published equations, parameter identity, calibration and solver accuracy remain distinct. A candidate calibration is usable for its explicitly defined model and can generate EOS while original-fit equivalence remains unresolved. An alternative energy functional cannot silently inherit the name and validation of another model.

## Policy correction

Earlier research/model-coverage reports used original-table fidelity and remaining native production ports among reasons to leave the overall goal partial. Current generation status is assessed from the actual physics route and its numerical support. Original-table discrepancies stay on the comparison axis. A Fortran production port is a separate engineering obligation when an independently checked Python physics implementation already generates native-readable tables. No past comparison or failed execution is rewritten as a success.

The 2026-11-16 base CI #618 passed both test suites and native calculations, but failed the root checksum manifest for data/research_generator_20261116/MODEL_COVERAGE.json. This is a packaging defect, not evidence against CLD physics; the current publication refreshes and checks actual artifact bytes.

## Core representations

- GM1: published rounded nonlinear RMF parameters, common nucleon rest energy939MeV, source-compatible npe. The physical functional and composition are the implementation target.
- APR: published21-parameter low/high functional, npemu and bulk two-charge Gibbs; Maxwell remains an explicit alternative. The author-code comparison is independent implementation evidence. Bulk Gibbs omits finite-size mixed-phase structure by definition; source equality is separately investigated.
- KDE0v1/SLy4: published full asymmetric Skyrme functionals with explicit kinetic/rest-mass conventions, density rearrangement and beta equilibrium. High-density causal crossings restrict physical use; they do not disable generation.
- BBB2: published Paris+Urbana SNM/PNM energies, the paper's parabolic asymmetry and a declared PCHIP representation. This reconstructs the published energy samples and solves composition; it does not re-solve the original Brueckner kernel. Total-energy thermodynamic pressure and the source-like baryon-path pressure are distinguished. The latter numerically explains the previous reference discrepancy but author-code attribution remains inferred.

## Crust representations

- SPG: source GDFM bulk, spherical CLD surface/Coulomb functional and explicitly identified surface fit. Compare complete stationarity/KKT, gas-species insertion, charge/baryon conservation, uniform competitors and density refinement. A finite-radius cutoff cannot stand in for a physical transition.
- GPPVA: local nucleon/electron densities and spatial sigma/omega/rho/Coulomb fields, local DD rearrangement where applicable, source field nonlinearities, Wigner–Seitz boundary conditions, cell-size and geometry minimization. A uniform bulk model plus fitted surface term does not replace this finite-range TF functional.
- Outer crust: discrete nuclear Gibbs minimization with source mass identity, electron/lattice corrections and neutron-drip/transition conditions. The historical common217-row component extends beyond pure outer crust; reproducing a pre-drip mass sequence does not fill its post-drip part.

TF profile composition and total species fractions are meaningful physical outputs. Discrete A,Z and bound/free decompositions require an explicit definition and may be unavailable for diffuse or infinite rod/slab geometries. Their absence must not be disguised as zero nuclei; companion profile data and native presence masks carry the distinction.

## Completion record

Report generation availability, paper/parameter conformity, numerical support, output/readback and reference agreement separately for every original model. Whole-project completion requires executable physical routes and their actual validated support; it cannot be inferred from test counts or source-file extraction. Remaining implementation and calibration questions receive concrete next calculations, not an unconditional CompOSE equality gate.
