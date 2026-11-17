# GPPVA spatial Thomas–Fermi regeneration

This implementation evaluates spatial finite-range RMF cells for GPPVA_DD2,
GPPVA_DDME2, GPPVA_TW, GPPVA_FSU2, GPPVA_FSU2R, GPPVA_FSU2H,
GPPVA_NL3wrL55 and GPPVA_TM1e. It is separate from the spherical CLD interface
and from historical proxy tables. The model fields/couplings are documented in
[GPPVA_BULK_FIELDS.md](GPPVA_BULK_FIELDS.md); source derivation and independent
functional review are in [TF_THEORY.md](TF_THEORY.md).
A source-table discrepancy is diagnostic, not an unconditional acceptance gate.

## Functional and constraints

The zero-temperature npe Thomas–Fermi approximation neglects shell corrections
and pairing. The radial coordinate is x=r/R in a Wigner–Seitz cell with dimension
d=1,2,3 and normalized volume weights proportional to x^(d-1). Flat, cylindrical
and spherical geometries permit slab, rod/tube and droplet/bubble profiles.
The matter variables are local neutron, proton and electron number densities in
fm^-3. Meson fields phi, omega and rho are in MeV, with hbar*c explicit.
The rho vertex is g_rho*tau_3/2. All fields have zero normal derivatives at both
ends; the weighted mean Coulomb potential V=e*A0 is zero.

The off-shell energy density in MeV fm^-3 is

    epsilon = <epsilon_kin,n + epsilon_kin,p + epsilon_kin,e
                 + U(phi)/(hbar*c)^3
                 + g_omega*omega*n_B + g_rho*rho*(n_p-n_n)/2
                 - W(omega,rho)/(hbar*c)^3>
              + G_phi - G_omega - G_rho + E_C,

where G_a=<a*(-laplacian_x)*a>/(2*hbar*c*R^2) and
E_C=<V*(n_p-n_e)>/2. Scalar energy is positive and vector terms negative
off shell; after solving the fields this yields the usual positive vector
interaction energy. Density-dependent couplings are evaluated locally.
The derivative includes the rearrangement

    Sigma_R = g_omega' * omega * n_B
            + g_rho' * rho * (n_p-n_n)/2 - g_sigma' * phi * n_scalar.

The discretized energy is minimized at fixed volume-averaged baryon density and
zero global charge. Electrons obey sqrt(k_Fe^2*c^2+m_e^2*c^4)-V=constant locally:
uniform electrons are not substituted. The electron PDE and neutrality are
solved inside the nucleon minimization. The reduced proton derivative is
mu_p(r)+mu_e,global; this follows from the electron envelope theorem.

The cell radius is optimized, with exact reduced-functional derivative

    d epsilon / d log(R) = -2 G_phi +2 G_omega +2 G_rho +2 E_C.

Consequently the virial identity at a stationary cell is
G_phi-G_omega-G_rho=E_C. Pressure is n_bar*mu_B-epsilon only after chemical
and radius stationarity. Independent nearby-density reoptimization checks this
pressure; a small optimizer objective change alone is not accepted as equilibrium.
Local electron Fermi energy is checked against the muon threshold, because this
specific crust model is npe.

## Discretization and numerical tests

The finite-volume radial Laplacian conserves the volume integral and is
self-adjoint in the volume-weighted inner product. The constant Neumann
zero mode is removed only for the electrostatic gauge. Scalar and vector
field equations use analytic Jacobians; the electron Newton solve uses
a scalar neutrality initialization and a gauge projection. The reduced matter
optimizer uses nonnegative local densities and a linear baryon constraint.
The explored radius range is 5–180 fm. Reaching a density or radius search
bound is recorded and prevents a nonuniform stationary result.

Acceptance for a finite candidate requires chemical KKT residual below
0.02 MeV and radius residual below 0.0002 MeV per baryon, together with
field equations, constraints, and successful optimization. These are numerical
criteria, not evidence of a global minimum. Mesh changes and neighboring density
samples are retained independently; a stationary but coarse cell is not promoted
to a continuum-converged EOS.

Tests cover a manufactured Neumann Poisson solution with second-order mesh
convergence, Fermi integrals against quadrature, uniform reduction of all eight
models, local density/radius derivatives, screened electron neutrality and
chemical equilibrium, and pressure against independent density differentiation.
Adversarial failures are preserved: weakly weighted central electron variables
could satisfy an energy stopping rule while violating chemical equilibrium;
exact electron elimination fixes that issue. A near-empty cached electron
chemical potential during an outer line search needed neutrality reinitialization.
A scalar neutrality bracketing endpoint also required outward rounding protection.

## Composition and native output

The calculated composition is the full set of local n_n(r), n_p(r), n_e(r),
plus their volume integrals. These include nucleons in diffuse clusters.
They are not automatically the CompOSE free-particle fractions Yn and Yp.
Finite A and Z are not assigned to rods/slabs, and a discrete bound/free partition
is not invented for diffuse spherical profiles.

The native table therefore contains thermodynamics, Ye, the declared zero Ymu,
and phase identifiers; free/cluster composition masks are unavailable in
nonuniform rows. Uniform rows contain free nucleon fractions. The companion
composition.json and selected_states.json.gz carry the physical total fractions,
profiles, geometry, radius and conservation errors. The representation is
spatial_TF_species_densities. This is available physical composition in a
different representation, not a thermo-only result.
compose_eos.tf_composition.read_tf_composition validates the companion against
the native table and profile conservation before returning total fractions.

These standalone cells do not establish a seam to a different outer-crust mass
functional. In particular, a constant energy shift cannot repair species-dependent
rest-mass conventions or change beta equilibrium consistently. Source BSk22
matching and full unified-crust claims require their own evidence.

## Reproduction

From the repository root, build the native executable and use:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=python_frontend \
    COMPOSE_EOS_FORTRAN_ROOT=. \
    COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran \
    python3 scripts/build_gppva_tf.py OUTPUT --workers 6

The output path must be new. All observed geometry seeds, failed optimizations,
mesh comparisons, native readback checks and source SHA-256 identities are retained.
A second command, scripts/refine_gppva_tf.py SOURCE OUTPUT, refines observed
competing branches to N192 and dilute cells to N384. Its selection is explicitly
limited to observed geometries within 0.01 MeV/baryon of the coarse minimum plus
the uniform reference. A finer-grid competitor outside that coarse window is
not mathematically excluded.

The finite campaign results and their supported intervals are recorded separately
in data/research_completion_20261116/gppva_tf. No stellar calculation or joined
core–crust EOS follows from this module.

## Executed campaign result

All eight models have 23/23 selected stationary source-interval nodes from 0.002 to 0.1 fm^-3. The final tables contain 226 rows including lower-density diagnostic extensions. Uniform phases at the upper end are explicitly identified. No held density is interpolated into the table.

| Model | Final rows | Lower-extension held densities (fm^-3) | Last clustered / first uniform sampled n | Max source-interval mesh P change |
|---|---:|---|---|---:|
| GPPVA_DD2 | 28 | 0.0002 | 0.07 / 0.075 | 0.1922% |
| GPPVA_DDME2 | 27 | 0.0002, 0.00035 | 0.07 / 0.075 | 0.1729% |
| GPPVA_TW | 28 | 0.0015 | 0.075 / 0.08 | 0.1724% |
| GPPVA_FSU2 | 28 | 0.0002 | 0.05 / 0.055 | 0.1780% |
| GPPVA_FSU2R | 29 | none | 0.09 / 0.095 | 0.3385% |
| GPPVA_FSU2H | 29 | none | 0.095 / 0.1 | 0.1704% |
| GPPVA_NL3wrL55 | 29 | none | 0.08 / 0.085 | 0.1000% |
| GPPVA_TM1e | 28 | 0.00075 | 0.095 / 0.1 | 0.2101% |


All 226 final profiles were independently reevaluated. Native write/readback and companion/table conservation checks succeeded for all eight models. The retained lower-density TF failures are diagnostic extensions below the source-used 0.002 fm^-3 interface, not holes in the regenerated source-used interval. Detailed receipts are in data/research_completion_20261116/gppva_tf/final and the corresponding gppva_bulk review.
