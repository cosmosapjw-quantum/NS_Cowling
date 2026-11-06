# Integrated compose_eos + relativistic Cowling pipeline (v3)

This release uses the already-generated thermodynamically consistent EOS tables as its default input. Unlike v2, the default path does **not** build a second q-linear (P(\mu_B)) EOS table and does not reintegrate pressure or rewrite energy-density knots.

## Read first

1. `docs/DIRECT_EOS_PROVIDER_V3.md`: physical derivation, interpolation decision, rejected alternatives, explicit source-knot preservation limits and validation.
2. `docs/TECHNICAL_MANUAL.md`: code architecture, units, TOV and Cowling equations.
3. `docs/ADVERSARIAL_AUDIT.md`: v2/v3 audit history and bounded corrections.
4. `docs/FINAL_INTEGRATION_RESULTS.md`: legacy cohort versus new direct-provider results.

## Default EOS representation

The provider uses a monotone PCHIP for (g=\epsilon-P) versus (\ln P). Original (P_i,\epsilon_i) knots are retained; (n_B(P)) is integrated from the cold first law and (\mu_B=(\epsilon+P)/n_B). This guarantees a continuous causal (0<c_s^2\le1) barotrope when the raw knot table passes the independent first-law and monotonicity gates. Source (n_B,\mu_B) columns are never overwritten; their reconstructed-knot differences are reported as interpolation systematics.

The v2 q-linear adapter remains available in `data/eos_cowling_adapter_v1` and `DATA_ADAPTERS_V2.zip` for comparison only. Its saved 234-row mode spectrum is a **legacy q-linear result**, not a v3 spectrum.

## Quick start

Install the bundled generator and Python package:

    python -m pip install -e eos_generator/compose_eos_pkg
    python -m pip install -e .

Audit the direct EOS representation without constructing replacement EOS CSVs:

    python scripts/compose_cowling.py audit-direct --source-dir data/eos_original_normalized --output-dir data/eos_direct_provider_v3

Run pipeline tests:

    python -m pytest -q pipeline/tests

Verify the old saved cohort without rerunning modes:

    python scripts/compose_cowling.py verify-existing --release oscillation_release_v1

The Fortran provider is in `oscillation_source/solutions/src/eos/ns_table_eos.f90`. The EOS callback is shared with the existing TOV and Cowling kernels. `cell_solve` uses the original table unchanged when the controlled-perturbation amplitude is zero.

## Scope

This package is confined to nonrotating, spherical, perfect-fluid (\ell=2) polar modes in the relativistic Cowling approximation, with cold beta-equilibrated barotropic perturbations and a finite positive pressure surface cutoff. It does not calculate metric perturbations or gravitational-wave damping, frozen-composition response, superfluid modes, sharp-interface conversion modes, or the exact zero-pressure surface.

The old REVTeX manuscript and 234-row spectrum describe the v2 q-linear-adapter cohort. The v3 direct-provider representative calculation is reported separately; the remaining 233 v3 spectra have not been recomputed.
