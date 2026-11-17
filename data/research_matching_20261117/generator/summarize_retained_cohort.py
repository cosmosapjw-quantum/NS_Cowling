"""Read and hash retained generation evidence; never rerun a physical solver."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]
DEST = Path(__file__).resolve().parent
BASE = ROOT / "data/research_native_20261116/public_cohort"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def grid(path):
    lines = path.read_text().splitlines()
    values = [float(x) for x in lines[2:]]
    assert int(lines[1]) - int(lines[0]) + 1 == len(values)
    return values


rows = []
for path in sorted(BASE.glob("*/attempt_001/COHORT_RESULT.json")):
    cohort = json.loads(path.read_text())
    evidence = cohort["native_evidence"]
    generated = path.parent / "output/generated"
    mismatches = [name for name, expected in evidence["output_sha256"].items()
                  if not (generated / name).is_file() or sha(generated / name) != expected]
    assert not mismatches, (cohort["model"], mismatches)
    nb = grid(generated / "eos.nb")
    phases = sorted({int(line.split()[3]) for line in (generated / "eos.compo").read_text().splitlines()})
    row = dict(model=cohort["model"], family=cohort["request"]["family"],
               historical_generation_status=cohort["status"],
               current_recorded_output_hash_verification="PASS",
               checked_output_files=len(evidence["output_sha256"]),
               retained_cohort=str(path.relative_to(ROOT)), retained_cohort_sha256=sha(path),
               rows=evidence["rows"], n_min=nb[0], n_max=nb[-1], phase_ids=phases,
               native_readback_tables=evidence["native_frontend_readback"],
               historical_failed_consistency_checks=evidence["native_consistency"]["generated"]["overall"]["failed_checks"],
               support=json.loads((generated / "SUPPORT.json").read_text()),
               support_sha256=sha(generated / "SUPPORT.json"))
    row["native_availability"] = "IMPLEMENTATION_VERIFIED_RETAINED_PUBLIC_NATIVE_EXECUTION"
    if row["family"] == "core":
        row["source_fidelity"] = {
            "APR_APR": "PUBLISHED_21_PARAMETER_FUNCTIONAL_WITH_EXPLICIT_BULK_TWO_CHARGE_GIBBS_ASSUMPTION",
            "BBB_BHF2": "PUBLISHED_SNM_PNM_ENERGY_RECONSTRUCTION; ORIGINAL_SMOOTH_FIT_DERIVATIVES_UNRESOLVED",
            "GM_GM1": "PUBLISHED_ROUNDED_RMF_PARAMETERS_COMMON_939_MEV_NUCLEON_REFERENCE_NPE",
            "RG_KDE0v1": "PUBLISHED_FULL_ASYMMETRIC_SKYRME_FUNCTIONAL_WITH_DECLARED_PARAMETER_PRECISION",
            "RG_SLY4": "PUBLISHED_FULL_ASYMMETRIC_SKYRME_FUNCTIONAL_WITH_DECLARED_PARAMETER_PRECISION",
        }[row["model"]]
        row["numerical_support"] = "RETAINED_REQUESTED_ROWS; DEFAULT_INTERVAL_NOT_A_UNIFORM_PHYSICAL_ACCEPTANCE_CLAIM"
        causal = {"APR_APR": .8689531515, "RG_KDE0v1": 1.4109658932, "RG_SLY4": 1.2650056715}
        if row["model"] in causal:
            row["approximate_causal_boundary_fm3"] = causal[row["model"]]
        if row["model"] == "BBB_BHF2":
            row["source_support_note"] = "Historical PCHIP at .7<n<=.8 depends on extrapolated .9 node; separate source-only9 reconstruction retained"
    elif row["family"] == "spg":
        row["source_fidelity"] = "EXPLICIT_CANDIDATE_SPHERICAL_GDFM_CLD; ORIGINAL_SURFACE_CALIBRATION_AND_VARIATIONAL_EXPORT_EQUIVALENCE_UNRESOLVED"
        row["numerical_support"] = "17_PUBLIC_ROWS; DENSER_RETAINED_NATIVE_CAMPAIGN_AND_STATE_DERIVATIVES_SEPARATE; NO_CROSS_SEGMENT_APPROVAL"
        inferred_path = DEST / "INFERRED_SPG_READBACK.json"
        inferred = json.loads(inferred_path.read_text())
        assert sha(ROOT / inferred["run_receipt"]) == inferred["run_receipt_sha256"]
        assert sha(ROOT / inferred["independent_summary"]) == inferred["independent_summary_sha256"]
        current = inferred["models"][row["model"]]
        row["conditional_inferred_variant"] = {"variant": "public_composition_inferred_20261116",
            "input_classification": "CONDITIONAL_PUBLIC_COMPOSITION_INFERENCE_NOT_ORIGINAL_AUTHOR_CALIBRATION",
            "new_native_generation_and_readback": current,
            "readback_receipt": str(inferred_path.relative_to(ROOT)),
            "readback_receipt_sha256": sha(inferred_path),
            "default_changed_by_this_lane": False}
    else:
        row["source_fidelity"] = "PUBLISHED_RMF_MODEL; INDEPENDENT_SPATIAL_TF_FINITE_SEARCH_NOT_COMPLETE_ORIGINAL_AUTHOR_CODE_EQUIVALENCE"
        row["numerical_support"] = "FINITE_DENSITY_GEOMETRY_SEED_MESH_SAMPLES; GLOBAL_MINIMUM_AND_MESH_CONTINUUM_UNPROVEN; CROSS_PHASE_INTERPOLATION_UNAPPROVED"
        row["nonuniform_composition_support"] = "ONLY_LEPTON_MASKS_PRESENT; FREE_CLUSTER_DECOMPOSITION_UNAVAILABLE"
    lower = generated / "lower_source_fit/eos.nb"
    if lower.exists():
        last = grid(lower)[-1]
        row["lower_last_density_fm3"] = last
        row["lower_to_inner_delta_n_fm3"] = nb[0] - last
        row["lower_to_inner_relation"] = "gap" if nb[0] > last else "overlap" if nb[0] < last else "same_density"
    dense = ROOT / "data/research_native_20261116/source_grid_cohort" / cohort["model"] / "attempt_001/COHORT_RESULT.json"
    if dense.exists():
        dense_cohort = json.loads(dense.read_text())
        dense_evidence = dense_cohort["native_evidence"]
        dense_dir = dense.parent / "output/generated"
        dense_mismatches = [name for name, expected in dense_evidence["output_sha256"].items()
                           if not (dense_dir / name).is_file() or sha(dense_dir / name) != expected]
        assert not dense_mismatches, (cohort["model"], dense_mismatches)
        row["retained_source_grid"] = dict(rows=dense_evidence["rows"],
            path=str(dense.relative_to(ROOT)), sha256=sha(dense),
            checked_output_files=len(dense_evidence["output_sha256"]), hash_verification="PASS")
    rows.append(row)
assert len(rows) == 18
result = dict(baseline="478c7b5de1514ac512f91156f36551478f36f688", purpose="REUSE_RETAINED_NATIVE_OUTPUTS_NO_PHYSICS_RERUN",
              claims={"native_generation_executed_now":False,"retained_output_hashes_checked_now":True,
                      "original_author_equivalence_complete":False,"continuous_support_certified_by_this_matrix":False},
              models=rows,
              totals=dict(models=len(rows),main_rows=sum(r["rows"] for r in rows),
                          native_tables=sum(len(r["native_readback_tables"]) for r in rows),
                          output_files_checked=sum(r["checked_output_files"] for r in rows),
                          source_grid_models=sum("retained_source_grid" in r for r in rows),
                          source_grid_rows=sum(r.get("retained_source_grid",{}).get("rows",0) for r in rows)))
(DEST / "RETAINED_COHORT_MATRIX.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result["totals"], sort_keys=True))
