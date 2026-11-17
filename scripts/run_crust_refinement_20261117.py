#!/usr/bin/env python3
"""Bounded native crust continuation, mesh sensitivity and actual sharp matching."""
from pathlib import Path
import argparse, csv, json, os, sys, time
from collections import Counter
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts"),str(ROOT/"pipeline")]
from compose_eos.tf_candidate import (run_state, interpolate_initial, export_compose, sha, save)
from compose_eos.api import diagnose_maxwell_phase, read_compose, write_compose
from compose_eos.maxwell_v2 import construct_maxwell_transition_v2, _segment_eosdata
from compose_eos.native import run_native, executable
from maxwell_three_family_campaign import read_prepared, valid_components, plain
from match_report_admission import admit_match_report
from match_core_crust import verify_internal_readback, verify_internal_metadata_readback
BASE=ROOT/"data/crust_refinement_20261117/tf_regeneration"
EXE=ROOT/"build/tf_candidate_G_final/tf_candidate_driver"
CONFIG={
 "DDME2":{"densities":np.linspace(.006,.014,9).tolist(),"seed_model":"DDME2","core":"new_native:PK1_CANONICAL:fine","mn":939.,"mp":939.},
 "TM1_CANONICAL":{"densities":np.linspace(.036,.052,9).tolist(),"seed_model":"TM1e","core":"native_generated:GM_GM1","mn":938.,"mp":938.}}
F_MANIFEST=ROOT/"data/new_model_maxwell_20261117/prepared/MANIFEST.json"

def load(path):
    return json.loads(Path(path).read_text())

def old_seed(model,density,dimension=3,seed="droplet"):
    folder=ROOT/f"data/research_native_20261116/source_grid_cohort/GPPVA_{model}/attempt_001/output/generated"
    selected=list(csv.DictReader((folder/"selected.csv").open()))
    node=min(selected,key=lambda s:abs(float(s["n_b"])-density))
    attempts=list(csv.DictReader((folder/"attempts.csv").open()))
    choices=[s for s in attempts if int(s["node"])==int(node["node"]) and int(s["dimension"])==dimension and s["seed"]==seed]
    good=[s for s in choices if s["status"] in ("STATIONARY","UNIFORM_STATIONARY")]
    chosen=min(good or choices,key=lambda s:float(s["epsilon"]))
    profile_rows=[s for s in csv.DictReader((folder/"profiles.csv").open())
         if int(s["node"])==int(node["node"]) and int(s["attempt"])==int(chosen["attempt"]) and int(s["mesh"])==int(chosen["mesh"])]
    profile_rows.sort(key=lambda s:int(s["cell"]))
    profile=[[float(s[k]) for k in ["nn","np","ne","mu_n","mu_p","mu_e","sigma","omega","rho","V"]] for s in profile_rows]
    assert len(profile)==int(chosen["mesh"])
    return {"id":f"legacy_{model}_node{node['node']}_attempt{chosen['attempt']}","n_b":float(node["n_b"]),
      "points":int(chosen["mesh"]),"dimension":dimension,"seed":seed,"radius":float(chosen["radius"]),
      "profile":profile,"legacy_status":chosen["status"],"epsilon":float(chosen["epsilon"]),
      "source_hashes":{str((folder/n).relative_to(ROOT)):sha(folder/n) for n in ["selected.csv","attempts.csv","profiles.csv"]}}

def solve(parent,model,n,points,label,dimension=None,seed=None):
    destination=BASE/"states"/model/label
    if (destination/"STATE.json").exists():
        s=load(destination/"STATE.json")
        assert s["model"]==model and s["points"]==points and abs(s["n_b"]-n)<1e-15 and s["executable_sha256"]==sha(EXE)
    else:
        nn,pp,r=interpolate_initial(parent,points,n)
        s=run_state(EXE,model,n,points,dimension or parent["dimension"],seed or parent["seed"],r,
           work_dir=destination,initial_nn=nn,initial_np=pp,parent=parent.get("state_path",parent.get("id")))
    # One fresh SLSQP initialization from the actually returned profile only.
    if not s.get("accepted_stationary") and "profile" in s and int(s.get("optimizer_mode",-1)) in (8,9):
        dest=destination.with_name(destination.name+"_warm")
        if (dest/"STATE.json").exists(): return load(dest/"STATE.json")
        nn,pp,r=interpolate_initial(s,points,n)
        warm=run_state(EXE,model,n,points,s["dimension"],s["seed"],r,work_dir=dest,
            initial_nn=nn,initial_np=pp,parent=s["state_path"])
        warm["productive_repair"]="one returned-profile optimizer-state restart; same ftol and gates"
        save(dest/"STATE.json",warm)
        s=warm
    print(json.dumps({k:s.get(k) for k in ["model","n_b","points","status","geometry","pressure","radius","wall_seconds"]}),flush=True)
    return s

def build_mesh(model,points,previous=None):
    cfg=CONFIG[model]; dens=cfg["densities"]; center=len(dens)//2
    rows=[None]*len(dens)
    seed=old_seed(cfg["seed_model"],dens[center]) if previous is None else previous[center]
    rows[center]=solve(seed,model,dens[center],points,f"N{points:03d}_i{center:02d}")
    for order in [range(center-1,-1,-1),range(center+1,len(dens))]:
        last=rows[center]
        for i in order:
            rows[i]=solve(last,model,dens[i],points,f"N{points:03d}_i{i:02d}")
            if "profile" in rows[i]: last=rows[i]
    save(BASE/f"{model}_N{points:03d}_STATES.json",rows)
    return rows

def connected_blocks(rows):
    blocks=[]; block=[]
    for s in rows:
        if s.get("accepted_stationary"):
            block.append(s)
        else:
            if len(block)>=3:blocks.append(block)
            block=[]
    if len(block)>=3:blocks.append(block)
    return blocks

def matched_component(low,high,lc,hc,dest):
    dest.mkdir(parents=True,exist_ok=False)
    out={"crust_component":lc,"core_component":hc,"status":"STARTED"}
    try:
        a=_segment_eosdata(low,(lc["phase_id"],lc["start"],lc["stop"]))
        b=_segment_eosdata(high,(hc["phase_id"],hc["start_zero_based"],hc["stop_exclusive"]))
        a.beta_equilibrated=True;b.beta_equilibrated=True
        matched,rep=construct_maxwell_transition_v2(a,b,assume_local_charge_neutrality=True,return_report=True)
        out["native_report"]=plain(rep.native_report);out["v2_report"]=plain(rep)
        matched.write_internal_directory(dest/"expected_internal")
        write_compose(dest/"compose",matched);back=read_compose(dest/"compose")
        back.write_internal_directory(dest/"readback_internal")
        out["readback"]=verify_internal_readback(dest/"expected_internal/eos_data.csv",dest/"readback_internal/eos_data.csv",expected_n=matched.n,actual_n=back.n)
        out["metadata"]=verify_internal_metadata_readback(dest/"expected_internal/eos_data.nml",dest/"readback_internal/eos_data.nml")
        out["admission"]=admit_match_report("maxwell",rep.native_report,dest/"readback_internal/eos_data.csv",dest/"readback_internal/eos_data.nml")
        out["consistency"]=run_native("verify",dest/"compose",dest/"CONSISTENCY.json")
        assert out["readback"]["status"]=="PASS" and out["metadata"]["status"]=="PASS"
        assert out["admission"]["status"]=="ADMITTED" and out["consistency"].get("overall")=="PASS"
        out["status"]="GENERATED_VERIFIED_FINITE_DOMAIN";out["compose_path"]=str(dest/"compose")
        out["table_sha256"]={p.name:sha(p) for p in (dest/"compose").glob("eos.*")}
    except Exception as exc:
        out.update(status="REJECTED",exception_type=type(exc).__name__,error=str(exc))
        if hasattr(exc,"report"):out["failure_report"]=plain(exc.report)
    save(dest/"RESULT.json",out)
    return out

def export_and_match(model,points,rows):
    cfg=CONFIG[model]; manifest=load(F_MANIFEST)
    core=next(s for s in manifest["source_records"] if s["id"]==cfg["core"])
    high=read_prepared(core)
    root=BASE/"tables_v2"/model/f"N{points:03d}"
    if (root/"SUMMARY.json").exists():return load(root/"SUMMARY.json")
    root.mkdir(parents=True,exist_ok=True)
    summary={"model":model,"points":points,"requested_densities":cfg["densities"],"accepted_requested_nodes":sum(s["accepted_stationary"] for s in rows),
       "core_source_id":cfg["core"],"core_prepared_manifest_sha256":sha(F_MANIFEST),"native_match_executable_sha256":sha(executable()),
       "scope":"finite observed branch; partner-independent fixed density window and native interval diagnosis; no global minimum claim","blocks":[]}
    for i,block in enumerate(connected_blocks(rows)):
        data,receipt=export_compose(root/f"block_{i:02d}",block,mn=cfg["mn"],mp=cfg["mp"],model=model)
        support={"scope":"computed_native_samples","global_minimum_proven":False,"continuous_mesh_limit_proven":False,
                 "cross_phase_interpolation":"UNAPPROVED","cross_segment_interpolation":"UNAPPROVED",
                 "physical_acceptance":"native stationarity and independent thermodynamic interval diagnosis"}
        save(root/f"block_{i:02d}/compose/SUPPORT.json",support)
        data.source_dir=str(root/f"block_{i:02d}/compose")
        diagnostic=diagnose_maxwell_phase(data);save(root/f"block_{i:02d}/NATIVE_INTERVAL_DIAGNOSTIC.json",diagnostic)
        record={"generation":receipt,"interval_rejections":dict(Counter(x["category"] for x in diagnostic["intervals"] if not x["valid"])),"components":[],"matches":[]}
        for j,(a,b) in enumerate(valid_components(diagnostic["intervals"],data.n)):
            lc={"start":a,"stop":b,"phase_id":int(data.Iphase[a]),"density":[float(data.nb[a]),float(data.nb[b-1])],
                "mu_E":[float((data.eps[a]+data.p[a])/data.nb[a]),float((data.eps[b-1]+data.p[b-1])/data.nb[b-1])]}
            record["components"].append(lc)
            for k,hc in enumerate(core["components"]):
                record["matches"].append(matched_component(data,high,lc,hc,root/f"block_{i:02d}/match_{j:02d}_{k:02d}"))
        summary["blocks"].append(record)
    save(root/"SUMMARY.json",summary)
    print(json.dumps({"matched":model,"mesh":points,"statuses":[m["status"] for b in summary["blocks"] for m in b["matches"]]}),flush=True)
    return summary

def main():
    ap=argparse.ArgumentParser();ap.add_argument("phase",choices=["production","refine128","derivative","reproject"])
    ap.add_argument("--models",nargs="+",default=list(CONFIG));args=ap.parse_args()
    BASE.mkdir(parents=True,exist_ok=True)
    if args.phase=="production":
        plan={"base_commit":"d89d877d50b6ae3ad4d5cdea29d5ff67e3a7d09a","models":CONFIG,
              "initial_meshes":[16,32,64],"warm_restart_limit":1,"fixed_ftol":1e-14,"fixed_maxiter":1800,
              "driver_sha256":sha(EXE),"registry_sha256":sha(ROOT/"fortran_eos/src/compose_eos_tf_models.f90"),
              "matching_binary_sha256":sha(executable()),"script_sha256":sha(Path(__file__)),
              "invariant":"fresh native physical states; no matching-target energy shifts or parameter fitting"}
        save(BASE/"INPUT.json",plan)
        for model in args.models:
            previous=None
            for points in [16,32,64]:
                rows=build_mesh(model,points,previous);export_and_match(model,points,rows);previous=rows
            anchor=CONFIG[model]["densities"][4]
            alternative=old_seed(CONFIG[model]["seed_model"],anchor,2,"droplet")
            alt=solve(alternative,model,anchor,64,"N064_anchor_rod",dimension=2,seed="droplet")
            save(BASE/f"{model}_N064_ROD_CONTROL.json",alt)
    elif args.phase=="reproject":
        for model in args.models:
            for points in [16,32,64]:
                export_and_match(model,points,load(BASE/f"{model}_N{points:03d}_STATES.json"))
    elif args.phase=="refine128":
        for model in args.models:
            rows=build_mesh(model,128,load(BASE/f"{model}_N064_STATES.json"));export_and_match(model,128,rows)
    else:
        for model in args.models:
            points=128 if (BASE/f"{model}_N128_STATES.json").exists() else 64
            rows=load(BASE/f"{model}_N{points:03d}_STATES.json");s=rows[4];n=s["n_b"];h=1e-5
            st={}
            for k in [-2,-1,1,2]:st[k]=solve(s,model,n+k*h,points,f"N{points:03d}_derivative_{k:+d}")
            checks={"model":model,"mesh":points,"n_b":n,"h":h,"center_state":s["state_path"],"states":{str(k):v["state_path"] for k,v in st.items()},
                    "all_accepted":all(v["accepted_stationary"] for v in st.values()),"same_geometry":all(v["geometry"]==s["geometry"] for v in st.values())}
            if checks["all_accepted"] and checks["same_geometry"]:
                ph=n*n*((st[1]["epsilon"]/st[1]["n_b"])-(st[-1]["epsilon"]/st[-1]["n_b"]))/(2*h)
                p2h=n*n*((st[2]["epsilon"]/st[2]["n_b"])-(st[-2]["epsilon"]/st[-2]["n_b"]))/(4*h)
                checks.update(euler_pressure=s["pressure"],P_h=ph,P_2h=p2h,
                    abs_Euler_difference=abs(ph-s["pressure"]),relative_Euler_difference=abs(ph-s["pressure"])/abs(s["pressure"]),
                    derivative_step_difference=abs(ph-p2h))
            save(BASE/f"{model}_THERMODYNAMIC_CHECK.json",checks)
    print("COMPLETED",args.phase,flush=True)
if __name__=="__main__":main()
