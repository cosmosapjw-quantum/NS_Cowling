# EOS Generator Scope and Explicit Limitations (compose_eos 3.3.1)

## Purpose

The `compose_eos` Python package is a cold, one-dimensional CompOSE EOS reader, thermodynamic auditor, table reconstructor, crust-core matcher, phase-transition utility and TOV helper. It is not, by itself, a universal microscopic EOS generator.

A successful parse, reconstruction, match, root search or serialization check does not establish that the underlying nuclear model is physically valid. Source admission and physical validation remain separate requirements.

## Supported input/output contract

The current `EOSData` and CompOSE I/O implementation support:

- One selected one-dimensional slice in baryon density (n_B).
- Exactly one temperature and one charge-fraction grid value per input table.
- Cold tables only: (T=0).
- CompOSE lepton flag (I_l\in\{0,1\}).
- Optional `Nadd=0) or `Nadd=1).
- Four particle IDs supported by the current data model.
- At most one averaged nuclear group with the supported group identifier.
- Internal nuclear-group baryon fraction (Y_h=A_{\rm av}Y_I), explicitly converted to/from the CompOSE group number abundance.
- Optional composition absence, kept distinct from known zero abundance.
- Model-specific phase identifiers retained when the corresponding composition file is available.

The current package does not implement a general multi-(T), multi-(Y_q), multi-species, arbitrary-(N_{\rm add}), multiple-nuclear-group CompOSE tensor reader/writer. Such inputs must be rejected or handled by a separate validated tool.

## Reconstruction limitations

`fill_density_gap` and `reconstruct_master` are table reconstruction/normalization procedures. They may regularize or project supplied fields under an explicitly selected prescription.

They do not independently regenerate the underlying RMF/EDF parameterized nuclear microphysics, nuclear masses, full cluster ensemble, surface tension, Coulomb lattice, or spatial density profiles.

A rebuilt table must retain a receipt describing all changed fields, source hashes, reconstruction method and residuals. Do not describe a projection or normalization as an independent microscopic EOS calculation.

## Crust microphysics

A complete compressible liquid-drop (CLD) or Thomas-Fermi (TF) crust generator generally requires a model-specific energy functional and a minimization problem including bulk, gradient/surface, Coulomb and lepton contributions, together with cell geometry and equilibrium constraints.

The current Python `compose_eos` package does not solve a general Wigner-Seitz TF spatial profile or a full pasta variational problem.

The repository has separate Fortran CLD/Wigner-Seitz components and tests. Their existence or successful synthetic unit tests do not establish an end-to-end production path. A production claim requires a native CLI/API route, real source input, deterministic output, independent read-back and per-case receipts.

Algebraic recovery of free neutron/proton fractions from a supplied averaged nuclear group is not equivalent to a spatial TF solution or variational minimization.

## Maxwell and Gibbs scope

### Maxwell

The sharp Maxwell path is for two cold, locally charge-neutral, beta-equilibrated one-dimensional phase barotropes. The input admission contract establishes the selected temperature and charge fraction, lepton inclusion, beta-equilibrium status, phase segmentation, thermodynamic consistency, stability and causality.

The Maxwell density and energy-density jumps must satisfy common pressure/common baryon chemical potential coexistence and the latent-energy identity.

### Gibbs

The Gibbs API requires two charge-resolved pressure surfaces (P(\mu_B,\mu_Q)) plus a lepton sector. It is not valid to substitute two one-dimensional beta-equilibrated barotropes.

For fixed (\mu_B), the current root isolator uses the union of both tensor-product spline knot vectors. It solves the quadratic charge-density difference on each knot span to locate all stationary points, then isolates pressure-equality roots on the resulting monotone intervals. This is complete for the declared piecewise-cubic interpolants; it does not recover arbitrary physics that is absent from the source grids.

The generated mixed-phase EOS is a thermodynamic barotrope. Species-resolved mixed-phase composition is not emitted by the current package. Global charge neutrality and volume fraction are recorded in the sidecar.

The surface inputs are supplied to the current API; the API does not itself generate those surfaces from an RMF, chiral EFT, quark-matter, or nuclear density-functional solver.

### CSS

Phenomenological constant-sound-speed branches are parameterized high-density models. They are not microscopic quark-matter calculations. A CSS bridge, a sharp Maxwell transition, and a smooth core-crust interpolation are distinct constructions and must retain distinct method identifiers and validation claims.

## Thermodynamic validation requirements

For a cold one-component beta-equilibrated barotrope:

[
d\epsilon=\mu_B\,dn_B,\qquad
P=n_B\mu_B-\epsilon,\qquad
\frac{dP}{d\mu_B}=n_B,
\]

[
c_s^2=\frac{dP}{d\epsilon}
=\frac{n_B}{\mu_B\,dn_B/d\mu_B}.
\]

Every interpolation or reconstruction must state which variable is authoritative, quantify deviations from source knots and validate positivity, mechanical stability and causality over intervals, not only at source rows.

## Remaining production extensions

The following capabilities remain separate extension work unless a specific native end-to-end receipt proves otherwise:

1. Multi-dimensional (T\times n_B\times Y_q) CompOSE tensor I/O.
2. Multiple nuclear groups and arbitrary composition species with row-level presence/applicability masks.
3. Provider-specific source semantics and model provenance.
4. Native Fortran CLD/TF generation from raw model parameters through final CompOSE output.
5. Full Wigner-Seitz spatial profiles and variational pasta minimization.
6. A unified RMF/EDF/CLD/TF solver-to-CompOSE generator interface.
7. Source-backed Gibbs surfaces and species-resolved mixed-phase composition.
8. Phase-interface-aware TOV and Cowling perturbation support.
9. Automated source-to-output regeneration with exact model, code and executable hashes.

Until those capabilities are implemented and validated, describe the current package as a bounded cold 1D EOS processing/matching tool with selected phase-transition APIs, not a universal multi-dimensional EOS generator.
