# Neutron-star quadrupolar f and p1 spectra (Cowling)

## Contents
- data/spectrum_summary_234.csv: one row per admitted joined EOS, with f and p1 frequencies, mass/radius, node counts, matching residuals and table hashes.
- data/eos_original_normalized/: 234 normalized source EOS tables, unchanged by the projection step.
- data/eos_projected/: 234 thermodynamically reconstructed EOS tables plus per-table projection receipts.
- data/representative_profiles/: nominal-tolerance background and matched inner/outer mode profiles for one representative EOS from each core family.
- paper/: REVTeX source, bibliography and figures; compiled PDF is at the release root.
- provenance/: source classifier, preparation manifest, matching metadata, native-cell audit, convergence evidence and source package archive.
- workflow/: compiled cell_solve, build manifest and scripts used in intake, projection, scan and assembly.

## Physical scope
All frequencies are ell=2 polar-fluid modes in the relativistic Cowling approximation, with a cold beta-equilibrated barotropic perturbation closure. The stellar surface is the lowest positive pressure in the EOS table. These results are not full-GR quasinormal-mode frequencies and do not include gravitational-wave damping, rotation, magnetic fields, elastic crust response, superfluid relative motion, frozen-composition perturbations, or sharp phase-interface modes.

## EOS projection
Raw source tables are preserved. The derived EOS representation uses monotone q(mu_B)=n_B/mu_B, n_B=mu_B q, P=P0+integral n_B dmu_B, and epsilon=mu_B n_B-P. This is an explicit thermodynamic interpolation reconstruction, not a silent modification of the source tables. See the per-EOS projection receipts and provenance/release_summary.json.

## Core-family medians and ranges (kHz)
- BBB_BHF2: f 2.5051 [2.4883, 2.5136], p1 7.0083 [6.6076, 7.1823]
- RG_SLY4: f 2.3698 [2.3574, 2.3803], p1 6.7177 [6.3471, 6.8906]
- GM_GM1: f 1.9202 [1.8893, 1.9401], p1 6.0090 [5.6866, 6.0773]

The frequency spread is a property of this 234-model cohort and should not be generalized to arbitrary neutron-star EOS families.
