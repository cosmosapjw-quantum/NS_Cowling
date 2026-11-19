#!/usr/bin/env python3
"""Compare native sharp Cowling break-point stopping against baseline.

Hold the EXACT same physical EOS/center frequency brackets while changing
only the numerical source waypoint policy and input tabulation resolution.
Do not call a method improvement without actual convergence evidence.
"""
from __future__ import annotations
import hashlib,json,sys,time,traceback
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import (NativeSharpCowlingDVODE,solve_sharp_background_native)
from compose_eos.stellar_cowling import CowlingSettings

BASE=ROOT/"data/cowling_maxwell_native_20261119"
OUT=ROOT/"data/cowling_eos_source_calibration_20261119/sharp_waypoints"
M03=json.loads((BASE/"validated_convergence/RESULT.json").read_text())
CASE_NAMES={
 "same_FSU2_maxwell":"three_layer_original_fsu2/maxwell/REPRESENTATION.json",
 "same_FSU2_css":"three_layer_original_fsu2/css/REPRESENTATION.json",
 "original_FSU2_css":"official_unified_css_sharp/FSU2/REPRESENTATION.json",
 "BBP_BPS_FSU2_css":"bbp_bps/whole_models/GPPVA_FSU2/css/REPRESENTATION.json",
}
SET=CowlingSettings(surface_h=3e-8,rtol=3e-9,atol=3e-12,max_step_h=.001)
TOV=dict(surface_h=3e-8,rtol=3e-10,atol=3e-11,max_step_h=.004)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+"\n")

def run_one(name,rel):
    p=BASE/rel
    source=M03["cases"][name]
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(p.read_text())))
    row=dict(name=name,source_path=str(p.relative_to(ROOT)),source_sha256=sha(p),
             original_4097_input_frequencies_hz=source["input_reference_frequencies_hz"],
             settings=dict(mode=vars(SET),tov=TOV),
             mode_results={})
    for resolution in (2049,4097,8193):
        start=time.monotonic()
        bg=solve_sharp_background_native(eos,.30,native_samples=resolution,**TOV)
        row["mode_results"][str(resolution)]=dict(background_M_msun=bg.mass_msun,
              background_R_km=bg.radius_km,source_ode_seconds=time.monotonic()-start,
              modes={})
        for break_policy in (False,True):
            solver=NativeSharpCowlingDVODE(bg,SET,samples=resolution,knot_waypoints=break_policy)
            tag="physical_break_waypoints" if break_policy else "baseline_continuous_DVODE"
            trials={}
            for conv in ("slow","rapid"):
                trials[conv]={}
                for mode in ("f","p1"):
                    f=float(source["input_reference_frequencies_hz"][conv][mode])
                    delta=.007*f
                    start=time.monotonic()
                    try:
                        z=solver.solve_mode((f-delta,f+delta),conversion=conv,
                                            profile_points=200)
                        if not (z["W_nodes_per_phase"]==z["V_nodes_per_phase"]==
                                {"f":0,"p1":1}[mode]):
                            raise RuntimeError("mode source node identity not retained")
                        trials[conv][mode]=dict(status="CHECKED",frequency_hz=float(z["frequency_hz"]),
                           root_residual=float(z["matching_residual"]),
                           join_relative=float(z["relative_join_mismatch"]),
                           source_break_waypoint_count=int(z["source_break_waypoint_count"]),
                           nfev=z["nfev"],elapsed_seconds=time.monotonic()-start)
                    except Exception as ex:
                        trials[conv][mode]=dict(status="REJECTED",error=str(ex)[:350],
                           elapsed_seconds=time.monotonic()-start)
            row["mode_results"][str(resolution)]["modes"][tag]=trials
            print("WAYPOINT",name,resolution,tag,
                  [(conv,mode,round(v.get("frequency_hz",-1),5),
                    round(v["elapsed_seconds"],3),v["status"])
                   for conv,rs in trials.items() for mode,v in rs.items()],flush=True)
        save(OUT/"PARTIAL.json",dict(last_completed=name,case=row))
    return row

def main():
    report=dict(schema="ns_cowling.sharp_knot_waypoint_convergence.v1",
       user_requested="controlled near-interface numerical accuracy/efficiency",
       interpretation="same physical EOS, TOV tolerance, frequency brackets and observed source content; vary only n_samples and optional forced eos-break stops",
       no_automatic_claim_of_improved_accuracy=True,
       cases={},status="RUNNING")
    for name,rel in CASE_NAMES.items():
        try:report["cases"][name]=run_one(name,rel)
        except Exception as ex:report["cases"][name]=dict(status="BLOCKED",type=type(ex).__name__,
                        error=str(ex),traceback=traceback.format_exc()[-1400:])
        save(OUT/"PARTIAL.json",report)
    report["status"]="REPRODUCIBLE_EXPERIMENT_COMPLETE"
    save(OUT/"RESULT.json",report)
    print("KNOT_CONVERGENCE_FINISH",[(k,"mode_results" in v) for k,v in report["cases"].items()],flush=True)
if __name__=="__main__":main()
