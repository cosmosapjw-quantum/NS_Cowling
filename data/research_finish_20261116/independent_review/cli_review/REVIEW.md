# Independent public generation CLI review

Base: `573a31485bc9af30915e01fe74f259e4fa5ac54f`. Read-only review of the new dispatcher and CLI, with isolated fixture mutations only.

Status: **PASS_SCOPED for inventory, request validation, retained core/SPG smoke, failure propagation, and the repaired TF evidence verifier.**

## Verified scope

- The catalog contains exactly 18 original identifiers: five core, five SPG and eight GPPVA TF models. Its region/representation contracts distinguish source functional implementations, the BBB published-energy reconstruction, candidate SPG surface calibration, and finite TF cell searches.
- Source-table equality is not an admission gate. Core requested intervals are preserved, while physical support remains a density-resolved diagnostic. SPG endpoints are independently solved; TF densities and radial resolutions must be explicit.
- Fresh output directories are required. Child errors and nonzero return codes survive in the dispatch receipt. Zero exit alone is insufficient; absent outputs, incomplete core grids, SPG HOLD and partial TF exports are rejected.
- The initial dispatcher test suite was actually rerun: **28 passed** (`existing_dispatch_tests.log`). The final expanded suite was independently rerun: **47 passed** (`final_dispatch_tests.log`).
- The retained actual CLI products were independently read with the native reader: GM1 3 rows, APR 5 rows from 3 requested densities plus two coexistence endpoints, and SPG M1 17 rows. Dimensional thermodynamics and species agreed with the retained scientific records to at most 1.14e-13. All declared scientific source hashes and the native executable SHA matched. See `CLI_SMOKE_REVIEW.json`.
- No new TF solver campaign was run by this reviewer.

## Reproduced defect

At dispatcher SHA `29dda811859f113785c94d0f2bafc693fc0ad6580e83a6bbe6e7c157e28ed3b1`, the TF verifier checked only nonempty native files, trusted export-status text, and profile/composition density lists. Unlike core/SPG, it did not enforce output hashes. A copied actual FSU2 N384 table was wrapped in truthful synthetic campaign metadata for its three rows. The verifier still accepted after its copied `eos.thermo` was replaced by plain text. It also accepted a selected row changed to `REQUIRES_REFINEMENT`.

The complete copied failure fixture, replay script, source SHA and observations are in `TF_EVIDENCE_RED.json`, `review_tf_evidence.py`, and `corrupted_tf_fixture/`. Only this isolated fixture was changed.

The final dispatcher binds trusted child readback to native/profile/composition bytes and interval metadata when present. It validates the retained selected-row stationary/feasibility/finite-thermodynamic contract, model and finest radial mesh, and every requested stage resolution. Independent replay accepted the clean copied actual FSU2 fixture and rejected all **12** byte-corruption, rehashed row-status/feasibility/model/mesh, missing-coverage and wrong-stage cases. See `TF_EVIDENCE_GREEN.json` and `review_tf_evidence_green.py`. Final dispatcher SHA: `e9b1634519be2795b3c77b3bac715eb1343920995f8ee489eecbcf01dee75267`; source bytes were unchanged during the replay. No new global-minimum or per-seed requirement was added.

The SPG documentation was corrected: its density grid ends at the clustered-side common-tangent endpoint, not at the separate fixed-density energy-crossing diagnostic.

## Limits

A dispatcher success certifies the requested finite output contract. It does not certify original author calibration, a complete crust, all-density physical support, continuous interpolation, spectral completeness, or a global energy minimum. The TF whole-run CLI path has no actual new solver smoke in this review; only its evidence contract is tested.

The green fixture supplies the required hash manifest explicitly; actual post-campaign builder manifest emission is owned by root integration and is not implied by this verifier-only replay. Two transient Remote Desktop Commander timeouts occurred while editing an independent replay script; a readback showed the edit had not executed, and it was then applied successfully. No production files were edited by this reviewer.
