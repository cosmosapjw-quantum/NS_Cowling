# Cowling phase-interface scope — sharp Maxwell backgrounds

**Repository:** `cosmosapjw-quantum/NS_Cowling`  
**Date:** 2026-11-15  
**Scope:** distinguish the general smooth-barotrope Cowling path from the separate phase-segmented Maxwell TOV/Cowling path, and determine what is still missing for EOS-generator core--crust CSS output.

## Background

At a cold first-order interface, pressure and baryon chemical potential are continuous while (n_B) and (epsilon) may jump, with
[
Deltaepsilon=mu_tDelta n_B.
]
Perturbation equations then require interface conditions. The physically correct conditions depend on the interface model and, for a fluid phase boundary, on whether phase conversion is rapid or slow.

## Implemented paths: do not conflate them

### General smooth-barotrope path

The general smooth-star solver under `oscillation_source/methods/` consumes a single smooth barotropic background. It does not automatically inspect a CompOSE table, discover an internal density jump, split the star into phases, or insert a phase-interface transfer.

A result obtained by passing a discontinuous table to this path is not thereby interface-validated.

### Separate sharp-Maxwell path

The repository also contains a dedicated phase-segmented implementation:

- `oscillation_source/solutions/src/structure/ns_maxwell_tov.f90` constructs a `maxwell_star_t` from separate high- and low-density EOS callbacks and validates coexistence at (P_t), (mu_B), and the cold first-law jump identity.
- `oscillation_source/solutions/src/modes/ns_maxwell_cowling.f90` integrates the core and envelope on separate pressure coordinates and applies an interface-transfer routine.
- `oscillation_source/solutions/src/modes/ns_interface_transfer.f90` provides distinct `material_interface_transfer` (slow conversion) and `rapid_interface_transfer` conditions.
- `oscillation_source/solutions/src/modes/ns_maxwell_mode_spectrum.f90` performs a phase-aware frequency scan and counts nodes separately on each side of the interface.
- `oscillation_source/solutions/tests_fortran/maxwell_checks.f90` and `maxwell_table_checks.f90` exercise analytic and table-backed Maxwell backgrounds, including slow/rapid Cowling shots and phase-aware mode-spectrum paths.

Therefore, the blanket statement that “the repository has no phase-interface Cowling code” is incorrect. The separate fluid-fluid Maxwell implementation exists and has native tests.

## What is still not established for core--crust CSS output

The EOS-generator route `scripts/match_core_crust.py` writes a thermo-only CompOSE table. It does not currently build the separate high/low EOS callbacks and `maxwell_star_t` required by `integrate_maxwell_tov`, nor does it run the resulting CSS crust-side jump through `shoot_maxwell_cowling`.

Consequently, there is no demonstrated end-to-end path:

[
	ext{matched CSS CompOSE output}
longrightarrow
	ext{phase-segmented TOV background}
longrightarrow
	ext{phase-interface Cowling spectrum}.
]

The existing Maxwell tests use analytic and normalized-table fluid-phase fixtures. They do not prove that a generated crust--CSS--core result is converted correctly into the separate Maxwell callbacks and passed through the full structure/mode path.

There is a second physical distinction: the current phase-segmented Maxwell solver treats both sides as barotropic fluids. It does not by itself implement an elastic solid crust with shear stresses or the solid-fluid crust-core boundary conditions. Thus it can support a fluid-fluid sharp interface under its stated conversion assumptions, but must not be advertised as a complete elastic crust-core oscillation treatment.

The single-speed CSS construction has a finite crust--CSS density jump and a continuous CSS--core join. Whether the crust can be represented by the fluid envelope callback is a model assumption that must be stated and tested, not silently inferred from the existence of the Maxwell solver.

## Status

- **IMPLEMENTATION_VERIFIED:** a separate sharp-Maxwell phase-segmented TOV/Cowling path exists; its native analytic/table fixtures test slow and rapid conversion.
- **UNRESOLVED:** automatic adapter from EOS-generator CSS output to the phase-segmented TOV/Cowling callbacks.
- **UNRESOLVED:** validation of the core--crust CSS interface under a realistic crust model, particularly if elastic crust physics is required.
- **NOT_CLAIMED:** full elastic crust perturbations or a universal crust-core interface prescription.

A smooth single-barotrope solver and the dedicated sharp-Maxwell solver are distinct numerical models. Do not flatten a density jump into the smooth solver's interpolation domain.

## Next bounded tasks

1. Build a tested adapter that splits a matched EOS at its exact one-sided interface rows and constructs separate low/high callbacks without interpolating across the jump.
2. Verify (P_t,mu_t,Delta n_B,Deltaepsilon), branch monotonicity, and first-law closure at the adapter boundary.
3. Run the generated background through the phase-segmented TOV routine and compare against independently assembled phase callbacks.
4. Run both slow- and rapid-conversion Cowling shots on that generated background.
5. Determine whether elastic crust physics is required; if so, derive and implement the solid-fluid boundary conditions rather than reusing fluid-fluid transfer formulas.
6. Add a real-data end-to-end fixture with source hashes, numerical residuals, and exact code revision.

## Literature scope

A sharp CSS density jump is a first-order transition parameterization, not a universal crust-core matching prescription. Crust-core matching methods differ in their physical assumptions and can affect neutron-star observables. A smooth unified crust-core EOS and a fluid-fluid Maxwell interface should be treated as distinct models.
