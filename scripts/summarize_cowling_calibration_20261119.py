#!/usr/bin/env python3
"""Collect source-locked COWLING-EOS-SOURCE-CALIBRATION-04 evidence.

Data report describes a scoped scientific result; it does not issue a
publication pass or convert model-dependent closures into original EOS.
"""
from __future__ import annotations
import hashlib,json,subprocess
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"data/cowling_eos_source_calibration_20261119"
M03=ROOT/"data/cowling_maxwell_native_20261119"
FILES={
 "independent_archive_bytes":"zenodo/ZENODO_INDEPENDENT_ARCHIVE_READBACK.json",
 "same_named_rmf":"source_audit/RESULT.json",
 "phase_local_thermodynamics":"phase_local_repairs/RESULT.json",
 "phase_local_stars":"phase_local_fixed_mass/RESULT.json",
 "coefficient_resolution":"resolution_separation/RESULT.json",
 "forced_eos_breakpoints":"sharp_waypoints/RESULT.json",
 "mode_family_overlap":"mode_continuation/RESULT.json",
 "tm1e_strict_conditional":"tm1e_mode_identity/RESULT.json",
 "tm1e_node_positions":"tm1e_mode_identity/NODE_LOCALIZATION.json",
}

def load(rel):return json.loads((D/rel).read_text())
def hashfile(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    receipts={key:dict(relative_path=path,sha256=hashfile(D/path))
              for key,path in FILES.items()}
    x={key:load(path) for key,path in FILES.items()}
    zen=x["independent_archive_bytes"]
    if not all(z["same_bytes_as_compose_download"] for z in zen.values()):
        raise RuntimeError("original ZIP remote independent archive byte identity failed")
    source={name:{
       "zenodo_record":zen[name]["record"],
       "observed_official_zip_sha256":zen[name]["observed_zip_sha256"],
       "official_listing_sha256_mismatch":zen[name]["compose_checksum_file_mismatch"],
       "provider_independent_bytes_verified":True,
       "core_same_density_comparison":{k:v for k,v in
            x["same_named_rmf"]["models"][name]["original_homogeneous_core_same_density_comparison"].items()
            if k not in ("samples",)},
       "original_unified_phase_audit":[{k:v for k,v in p.items() if k!="first_law"}
              |{"first_law_max_abs_F_shift_MeV":p["first_law"].get("max_abs_delta_epsilon_per_baryon_MeV"),
                "valid_on_full_phase":p["first_law"].get("valid_pressure_primary")}
              for p in x["same_named_rmf"]["models"][name]["phase_contiguous_audits"]]
      } for name in ("FSU2","TM1e","DDME2")}
    phase={}
    for name,item in x["phase_local_thermodynamics"]["cases"].items():
        p=item["variants"]["phasewise_firstlaw_integrated"]
        phase[name]=dict(maxwell_native_status=p["maxwell"]["status"],
               css_native_status=p["css"]["status"],
               transition_pressure_MeVfm3=p["css"].get("transition_pressure_MeV_fm3"),
               transition_density_fm3=[p["css"].get("n_low"),p["css"].get("n_high")],
               center_density_030_tov_status=p["css"].get("stellar",{}).get("status"))
    fixed={}
    for name,c in x["phase_local_stars"]["cases"].items():
        fixed[name]=dict(status=c["status"],
                 n_c_fm3=c.get("native_TOV",{}).get("nc_fm3"),
                 R_1p4_km=c.get("native_TOV",{}).get("R_km"),
                 delta_R_vs_cross_phase_repair_km=c.get(
                     "phasewise_minus_cross_seam_firstlaw",{}).get("R_km"),
                 delta_slow_p1_vs_cross_phase_repair_Hz=c.get(
                     "phasewise_minus_cross_seam_firstlaw",{}).get(
                         "mode_frequency_hz",{}).get("slow",{}).get("p1"),
                 conversion_mode_status={k:v["status"] for k,v in c.get(
                     "native_cowling",{}).items()})
    err=x["coefficient_resolution"]["cases"]["same_FSU2_maxwell"]
    original=json.loads((M03/"validated_convergence/RESULT.json").read_text())
    python_benchmark=float(original["cases"]["same_FSU2_maxwell"][
        "python_DOP853_oracle"]["modes"]["slow"]["p1"]["frequency_hz"])
    p4097=err["python_background"]["coefficient_sampling"]["4097"]["modes"]["slow"]["p1"]["f_Hz"]
    p32769=err["python_background"]["coefficient_sampling"]["32769"]["modes"]["slow"]["p1"]["f_Hz"]
    coarse_err=abs(p4097-python_benchmark)
    fine_err=abs(p32769-python_benchmark)
    if not (fine_err<coarse_err and fine_err<.005):raise ValueError("actual Cowling coefficient refinement failed")
    comparisons=[]
    times=[]
    for v in x["forced_eos_breakpoints"]["cases"].values():
        for resolution,trial in v["mode_results"].items():
            baseline=trial["modes"]["baseline_continuous_DVODE"]
            forced=trial["modes"]["physical_break_waypoints"]
            for conv in ("slow","rapid"):
                for label in ("f","p1"):
                    a=baseline[conv][label];b=forced[conv][label]
                    if not (a["status"]==b["status"]=="CHECKED"):
                        raise ValueError("unpaired original source waypoint method")
                    comparisons.append(abs(a["frequency_hz"]-b["frequency_hz"]))
                    times.append(b["elapsed_seconds"]/a["elapsed_seconds"])
    overlaps=[]
    for model in x["mode_family_overlap"]["cases"].values():
        if model["status"]!="SCOPED_MODE_FAMILY_CONTINUED":
            raise ValueError("model eigenshape continuation incomplete")
        for rows in model["conversion"].values():
            for row in rows:
                v=row.get("continuation_from_previous")
                if v:
                    if v["status"]!="CONTINUED":raise ValueError("bounded family discontinuity")
                    overlaps.extend(p["absolute_overlap"] for p in v["pairs"])
    tm1e=x["tm1e_strict_conditional"]["cases"]["slow"]
    jumps=[v for v in tm1e if v["continuation"] and
            v["continuation"]["node_pair_safe_tracker"]["status"]=="UNRESOLVED"]
    if not (len(jumps)==2 and all(v["continuation"][
        "conditional_shape_continuity_tracker"]["status"]==
        "CONTINUED_WITH_V_NODE_CHANGE" for v in jumps)):
        raise ValueError("conditional V node audit failed")
    out=dict(schema="ns_cowling.cowling_source_calibration_stage04.v1",
       date_utc_note="report generated on the authenticated computation workstation",
       parent_03_commit="0d802d09fff57df1cc5557916be824cbc457d2b7",
       current_head_at_evidence_capture=subprocess.check_output(
          ["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
       evidence_receipts=receipts,
       original_eos_provider_identity=source,
       phase_local_reconstruction_css=phase,
       same_mass_1p4_phase_local_stars=fixed,
       coupled_error_budget=dict(
           same_fsu2_maxwell_nc_fm3=.30,
           DOP853_Python_mode_reference_f_Hz=python_benchmark,
           python_TOV_native_Cowling_p1_4097_Hz=p4097,
           python_TOV_native_Cowling_p1_32769_Hz=p32769,
           baseline_abs_error_Hz=coarse_err,
           fine_abs_error_Hz=fine_err,
           approximate_reduction_factor=coarse_err/fine_err,
           principal_numeric_error_mechanism="native Cowling source EOS coefficient sampling, not native TOV"),
       forced_breakpoint_experiment=dict(
           comparisons=len(comparisons),
           max_absolute_frequency_change_Hz=max(comparisons),
           median_single_run_time_ratio_forced_over_baseline=float(np.median(times)),
           default_policy="off, no frequency accuracy benefit observed"),
       mode_family_continuation=dict(
           models=len(x["mode_family_overlap"]["cases"]),
           eigenfunction_transitions=len(overlaps),
           minimum_sign_independent_overlap=float(min(overlaps)),
           maximum_sign_independent_overlap=float(max(overlaps)),
           TM1e_strict_unresolved_V_node_changes=len(jumps),
           TM1e_conditional_shape_changes=len(jumps),
           source_both_side_zero_localization_support=True),
       BBP_primary_original_digital_source_status="NOT_VERIFIED_ONLY_Diener_2008_SECONDARY_35_ROW_TABLE",
       publication_gate="NOT_CLAIMED",
       microscopic_unification_gate="NOT_REQUIRED_FOR_EFFECTIVE_THREE_LAYER_TOV_COWLING",
       implementation_validation="PENDING_FRESH_PYTEST",
       status="SCOPED_CALIBRATION_RECEIPTS_GENERATED")
    # Bind reported implementation verification to source SHA exactly as
    # observed during the fresh test execution. Never reuse stale test logs.
    receipt=D/"verification/VERIFICATION_RECEIPT.json"
    if receipt.is_file():
        verified=json.loads(receipt.read_text())
        source_unchanged=all(hashfile(ROOT/rel)==digest for rel,digest in
          verified["exact_tested_source_sha256"].items())
        valid=(source_unchanged and verified["pytest"]["exit_status"]==0
               and all(x["exit_status"]==0 for x in
                       verified["Fortran_compile_with_bounds_checks"]))
        out["implementation_validation"]=dict(
            status="PASS_SOURCE_LOCKED_FRESH" if valid else
                   "STALE_OR_FAILED_SOURCE_VALIDATION",
            pass_count=verified["pytest"]["passed"],
            skipped_count=verified["pytest"]["skipped"],
            exact_tested_source_SHA_match=bool(source_unchanged),
            receipt_sha256=hashfile(receipt))
        if valid:
            out["status"]="SCOPED_CALIBRATION_SOURCE_AND_IMPLEMENTATION_VERIFIED"
    target=D/"CALIBRATION_OUTCOME.json"
    target.write_text(json.dumps(out,indent=2)+"\n")
    print("CALIBRATION_OUTCOME",target,"source_models",list(source),
          "phase_local",[(k,v["css_native_status"]) for k,v in phase.items()],
          "overlap_transitions",len(overlaps),"fractional_factor",coarse_err/fine_err,flush=True)

if __name__=="__main__":main()
