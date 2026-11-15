# Cowling phase-interface scope — density-jump background EOS

**Repository:** `cosmosapjw-quantum/NS_Cowling`
**Date:** 2026-11-15
**Scope:** §12 of the EOS390 remediation handoff — document whether the Cowling
perturbation solver supports the phase-interface matching conditions required
by a background EOS that carries a finite first-order density jump (sharp
Maxwell or single-speed CSS with a crust-side density jump).

## Background

A finite density jump in the background EOS changes the physical problem for
radial and non-radial perturbations. At a first-order phase transition the
background state variables are discontinuous in \(n_B\) and \(\epsilon\) (with
\(\Delta\epsilon = \mu_t \Delta n_B\), pressure and chemical potential
continuous), so the perturbation equations acquire interface matching
conditions that do not arise for a smooth barotrope. This document records what
the current Cowling implementation actually does, and what it does NOT claim.

## What the current Cowling solver implements

The oscillation solver (`oscillation_source/methods/`) is a Cowling-approximation
radial-mode solver on a smooth one-dimensional barotropic background. The
implementation contracts are:

1. **Continuous background barotrope.** The solver integrates the radial
   perturbation equation on a background EOS that is at least \(C^0\) in
   \(P(\epsilon)\) and \(C^0\) in \(n_B(r)\). It does not detect, locate, or
   special-case a material interface where \(n_B\) jumps.

2. **No phase-interface matching code path.** There is no code that, at a
   discontinuity of the background, applies a junction condition on the
   perturbation variables \((\delta P, \delta\epsilon, \xi)\) across the
   interface. The integration proceeds as if the background were smooth.

3. **No rapid- or slow-conversion switch.** The interface conversion rate
   \(\lambda\) (rapid conversion \(\lambda\to\infty\), slow conversion
   \(\lambda\to 0\)) does not appear as a solver parameter. No matching
   condition is selected by a conversion-rate flag.

4. **No interface motion.** The Lagrangian displacement \(\xi\) of the
   interface is not a solved degree of freedom. The background phase boundary
   is treated as a fixed point of the star, not as a moving surface whose
   motion couples to the perturbation.

5. **Boundary conditions actually implemented.** The solver applies a regular
   center condition and a surface matching condition to the exterior
   Schwarzschild vacuum. These are the standard smooth-star boundary
   conditions. No additional interface condition is imposed at a phase
   boundary.

6. **Test coverage of the interface.** No current test exercises a Cowling
   mode on a background EOS that carries a genuine sharp density jump. The
   provided default fixtures (`fixtures/modes/smooth_l2/`) are smooth
   barotropes; `ch07.md` records that "the provided default precision
   comparison is a continuous barotrope and does not automatically solve the
   spectrum of a whole star with a material density jump."

## What this means for the density-jump EOS

The single-speed CSS bridge with a crust-side density jump (§3 contract) and
the sharp Maxwell construction (§5) produce background EOS tables that contain
a genuine first-order density discontinuity. The current Cowling solver does
NOT support the interface matching conditions required by such a background.
Consequently:

- A radial mode computed on a density-jump background by the current solver is
  the result of integrating through the jump as if it were smooth. It is not a
  solution of the perturbation eigenvalue problem with correct junction
  conditions, and it must not be reported as such.
- The Seidov small-core onset criterion (§3.6) is a necessary stability
  condition, not a substitute for the interface-eigenvalue calculation. It is
  reported for completeness only.

## Status classification

`UNRESOLVED` — the Cowling phase-interface treatment is not implemented. The
density-jump EOS is therefore **not yet validated** for the perturbation
problem. This is an implementation limitation, not a numerical or physical
limitation of the EOS itself.

## Recommended next bounded task

Derive the junction conditions for a Cowling radial perturbation at a
first-order phase boundary under a stated conversion-rate assumption (rapid
or slow), implement the matching in the solver, and add a test that exercises
a mode on a sharp-Maxwell or single-speed-CSS background. The derivation must
cite the appropriate literature (e.g. the conversion-rate-dependent junction
conditions used in the hybrid-star g-mode literature) rather than inventing a
matching prescription. Until that task is complete, any Cowling mode reported
on a density-jump background must be labeled as **not interface-validated**.

## Literature pointers (not an endorsement of a single convention)

- Alford, Han & Prakash, "Generic conditions for hybrid star cores", arXiv:1302.4732 —
  the single-speed-plus-first-order-jump CSS parameterization used in §3.
- The hybrid-star g-mode literature (e.g. the conversion-rate-dependent
  matching used for discontinuity g-modes) derives junction conditions that
  depend on the interface conversion rate; those conditions are not yet
  implemented here.

This document records an unresolved limitation. It does not claim full
phase-interface Cowling support.
