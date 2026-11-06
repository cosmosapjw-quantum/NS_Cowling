# EOS generation to relativistic Cowling oscillations

> v3 authority note: The default provider is now documented in DIRECT_EOS_PROVIDER_V3.md. Sections below describing the q-linear P(mu_B) adapter are the historical v2 comparison path, not the v3 default. The source EOS tables remain authoritative and unmodified.

## 1. What is generated, and what the oscillation solver actually consumes

The upstream generator is the bundled compose_eos 3.3.0 release. Its thermodynamic reconstruction uses a master free-energy-per-baryon function F(n_B), with x=ln(n_B) and Q_1=dF/dx. At zero temperature and fixed beta-equilibrated composition it obtains

\[
P=n_B Q_1,\qquad \mu_B=F+Q_1,\qquad \epsilon=n_B F.
\]

Consequently, the knot values obey

\[
\epsilon+P=\mu_B n_B,\qquad d\epsilon=\mu_B\,dn_B,\qquad dP=n_B\,d\mu_B.
\]

The crust-core matching code likewise constructs smooth bridges from common thermodynamic variables. These statements concern the source reconstruction and the resulting table knots. They do not, by themselves, specify a unique continuous function between every pair of saved knots.

The supplied oscillation code makes a further, explicit representation choice. Its native interpolator is a cubic Hermite polynomial P(mu_B) on every table interval, with endpoint values P_i and slopes P'(mu_i)=n_i. Within that polynomial it defines n=P', epsilon=mu_B P'-P and chi=dP/d epsilon=P'/(mu_B P''). It then analytically checks each entire cell for positive compressibility and causality. This is not an EOS-generator check at the knots. It is an admissibility check on the particular continuous interpolation used by the differential-equation solver.

The raw joined tables retain their first-law identity at the knots, yet the raw-knot Hermite representation did not pass the chosen full-cell admissibility gate in the earlier analysis. The additional table is therefore a solver-compatible interpolation closure. It is not a new microscopic nuclear model, not a new crust-core matching prescription, and not a physical phase transition.

## 2. The additional q-linear representation

Set

\[
q(\mu_B)=\frac{n_B(\mu_B)}{\mu_B}.
\]

The construction takes q to be linear on each chemical-potential interval. Its endpoint values are the source q_i unless the optional weighted isotonic regression is needed to enforce nondecreasing q. Then

\[
n_B(\mu_B)=\mu_B q(\mu_B),
\]
\[
P(\mu_B)=P(\mu_{B,0})+\int_{\mu_{B,0}}^{\mu_B}xq(x)\,dx,
\]
\[
\epsilon(\mu_B)=\mu_B n_B(\mu_B)-P(\mu_B).
\]

Since q is linear in each interval, n_B is quadratic and P is exactly cubic. The result is therefore exactly represented by the native endpoint cubic-Hermite P(mu_B) interpolator, not merely approximated by it.

Differentiating gives

\[
\frac{dP}{d\mu_B}=n_B,\qquad
\frac{d\epsilon}{d\mu_B}=\mu_B\frac{dn_B}{d\mu_B}.
\]

Since n_B=mu_B q,

\[
c_s^2=\frac{dP}{d\epsilon}
=\frac{n_B}{\mu_B\,dn_B/d\mu_B}
=\frac{q}{q+\mu_B q'}.
\]

For positive mu_B and q, a nonnegative q' gives 0<c_s^2<=1. It also gives positive n'_B, epsilon' and P'. This is the precise reason for the q monotonicity constraint. No speed of sound is clipped after the fact.

The exact cell integral for h=mu_{i+1}-mu_i and delta_q=q_{i+1}-q_i is

\[
P_{i+1}-P_i
=h\left[\left(\mu_i+\frac h2\right)q_i+
\left(\frac{\mu_i}{2}+\frac h3\right)\delta_q\right].
\]

This is the expression used by the builder and is evaluated in double precision.

## 3. What happened in the 234-EOS cohort

The audit compared the prepared source knots with their q-linear solver adapters. All 234 source q sequences were already nondecreasing. Weighted PAV therefore pooled zero knots: it was a no-op to roundoff for this particular cohort. The effective operation was not a change in q or n_B. It was the choice of a piecewise-linear q between the existing knots, followed by pressure reintegration and the corresponding energy-density reconstruction.

Recorded cohort-wide changes are:
- maximum relative change in n_B: 3.1769e-16, roundoff level;
- maximum absolute pressure change divided by P_max: 1.2931e-3; median 5.0781e-4;
- maximum absolute energy change divided by mu_B n_B: 8.3590e-4; median 3.6598e-4;
- adapter first-law residual: at most 2.163e-16.

The source tables and adapter tables are retained separately. Source knots are never overwritten. The output receipt records both hashes, the projection blocks, pressure/energy differences, the first-law residual and the exact interpolation construction. For other EOS inputs where q is not monotone, PAV can change q and n_B; that is a real alteration of the tabulated thermodynamic curve and must be treated as a model change, not as a harmless formatting repair. The production command records this separately and fails if the resulting table is not positive, stable and causal.

## 4. Why first-law consistency does not imply interpolation admissibility

The first law constrains derivatives of one continuous thermodynamic curve. A finite table stores only samples of that curve and some derived quantities. An interpolator adds a further mathematical choice. A cubic Hermite interpolant fixed by P_i, n_i, P_{i+1}, n_{i+1} can overshoot between knots even when every endpoint satisfies epsilon_i+P_i=mu_i n_i. The generated curve can then have P''<=0 or mu_B P''-P'<0 somewhere inside a cell, which violates positive compressibility or c_s^2<=1 in the native representation.

This distinction is central:
1. source EOS construction and crust-core matching establish the tabulated thermodynamic values;
2. the interpolation contract defines the continuous function used by the solver;
3. the full-cell analytic gate decides whether that chosen function is mechanically stable and causal;
4. the q-linear adapter selects a different continuous function that satisfies the solver's P(mu_B) Hermite contract.

No claim is made that every possible interpolation of the raw table is inadmissible. The result is representation-specific. Alternative closures, such as a different monotone thermodynamic interpolant, are valid research options but must be compared as separate EOS realizations.

## 5. Provenance and reversibility

There are three non-interchangeable objects:
- source CompOSE and matched-EOS records from the EOS generator;
- normalized source-knot CSVs obtained from the successful matched tables;
- q-linear Cowling-adapter EOS tables, each accompanied by an adapter receipt.

The adapter must never replace the source CompOSE files. The source table hash, source matching metadata hash, adapter output hash, first-law residual and maximum changes are retained. Rebuilding the adapter with the bundled pipeline should reproduce the saved tables to floating-point serialization tolerance.

## 6. Reproduction

Install the bundled EOS generator and pipeline in the integrated release root:

    python -m pip install -e eos_generator/compose_eos_pkg
    python -m pip install -e .

To prepare source-knot tables from a previously generated matched-EOS directory:

    python scripts/compose_cowling.py prepare --matched-root /path/to/matched_eos_v3_2 --classifier eos_generator/docs/MATCHING_FAILURE_CAUSE_CLASSIFICATION_v3_3.json --output-dir data/eos_original_normalized

To create the q-linear adapter:

    python scripts/compose_cowling.py build-adapters --prepared-dir data/eos_original_normalized --output-dir data/eos_cowling_adapter

To verify the already completed f and p1 spectrum without running any integrations:

    python scripts/compose_cowling.py verify-existing --release oscillation_release_v1

The existing f/p1 results are included in this package. Verification is read-only. No mode calculation is started implicitly. The original Fortran solver, source equations, input namelists, legacy workflow scripts, results and validation records are included under oscillation_source and oscillation_release_v1.

## 7. Model scope

The oscillation calculation is the nonrotating, spherical, perfect-fluid, ell=2 polar problem in the McDermott-Van Horn-Scholl relativistic Cowling prescription: all Eulerian metric perturbations are neglected. The fluid perturbation uses the cold beta-equilibrated barotropic derivative chi=dP/d epsilon. It does not compute full-GR quasinormal modes, spacetime w modes or gravitational-wave damping times. Rotation, magnetic stresses, crust elasticity, superfluid relative motion and frozen-composition response are outside this model.

The finite positive pressure cutoff is part of the numerical model. It is not an exact vacuum surface. The matching label “maxwell” in the original 390-member matrix is a legacy smooth bridge label; it does not insert a sharp first-order density jump. Discontinuity g modes require a genuine sharp interface and conversion-rate-dependent junction conditions; they are not contained in these smooth one-dimensional barotropes.
