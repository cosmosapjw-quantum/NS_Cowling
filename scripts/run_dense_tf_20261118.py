"""Independent spatial-TF N128 subgrid states from frozen accepted H/G profiles.

Each attempt solves the unchanged native field equations at a NEW density
between H source knots. No interpolated estimate can become reference truth.
Restrict comparison to identical droplet geometry / branch / mesh.
"""
from pathlib import Path
import hashlib,json,sys,time,argparse,traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.tf_candidate import run_state,interpolate_initial,sha
EXE=Path("/home/cosmosapjw/ns_cowling_wt/scratch/native_dense_interp_20261118/tf_build/tf_candidate_driver")
OUT=ROOT/"data/native_dense_interpolation_20261118/tf"
CASES={
 "DDME2":dict(model="DDME2",n_low=.00725,n_high=.01075,step=.00025,
  src="data/crust_refinement_20261117/tf_regeneration/states/DDME2/N128*/**/STATE.json",
  branch="G/H DDME2 N128 droplet"),
 "DDMEX":dict(model="DDMEX_ORIGINAL_TABLE_RHO0",n_low=.007125,n_high=.009875,step=.00025,
  src="data/low_density_crust_20261117/low_extension/states/DDMEX_ORIGINAL_TABLE_RHO0/N128*/**/STATE.json",
  branch="H DDME-X original-rho0 N128 droplet"),
 "NL3":dict(model="NL3_PROVIDER_PRECISION",n_low=.02925,n_high=.03275,step=.0005,
  src="data/low_density_crust_20261117/low_extension/states/NL3_PROVIDER_PRECISION/N128*/**/STATE.json",
  branch="H provider NL3 N128 droplet"),
 "NL3wrL55":dict(model="GPPVA_NL3wrL55",n_low=.010125,n_high=.012875,step=.00025,
  src="data/low_density_crust_20261117/support_refinement/states/GPPVA_NL3wrL55/continued_G_lower_adjacent_branch/N128/**/STATE.json",
  branch="H NL3wrL55 lower adjacent branch N128 droplet"),
 "TM1e":dict(model="GPPVA_TM1e",n_low=.007875,n_high=.010125,step=.00025,
  src="data/low_density_crust_20261117/support_refinement/states/GPPVA_TM1e/continued_G_upper_F_component/N128/**/STATE.json",
  branch="H TM1e upper F branch N128 droplet; includes failed original knots")
}
def target_densities(spec):
    size=int(round((spec["n_high"]-spec["n_low"])/spec["step"]))+1
    return np.round(np.asarray(spec["n_low"])+np.arange(size)*spec["step"],9)
def existing_anchors(spec):
    rows=[]
    for path in sorted(ROOT.glob(spec["src"])):
        try:
            state=json.loads(path.read_text())
            if state.get("model")==spec["model"] and state.get("accepted_stationary") and state.get("points")==128 and state.get("geometry")=="droplet":
                state["_source_path"]=str(path.relative_to(ROOT))
                state["_source_sha256"]=sha(path)
                rows.append(state)
        except (OSError,ValueError,json.JSONDecodeError):
            continue
    return rows
def choose_anchor(states,density):
    result=min(states,key=lambda x:(abs(x["n_b"]-density),x["_source_path"]))
    if abs(result["n_b"]-density)>.00175:
        raise ValueError("no nearby source same-branch profile; refuse blind seed")
    return result
def save(path,x):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def compute(case):
    if case not in CASES:raise ValueError("case not supported")
    spec=CASES[case]
    target=target_densities(spec)
    source=existing_anchors(spec)
    if not source:raise RuntimeError("no source N128 accepted profile: "+case)
    dest=OUT/case
    if dest.exists():raise FileExistsError("existing dense TF run, preserve it: "+str(dest))
    dest.mkdir(parents=True)
    rec=dict(schema="ns_cowling.native_dense_tf_interleaved.v1",case=case,
        model=spec["model"],source_branch=spec["branch"],
        target_densities_fm3=target.tolist(),
        actual_native_executable_sha256=sha(EXE),
        physical_acceptance="Original TF run_state accepted_stationary, strict optimizer and constraints; geometry droplet; no gap filling",
        status="RUNNING",records=[])
    save(dest/"RESULT.json",rec)
    for i,n in enumerate(target):
        t=time.monotonic()
        print("TF_ATTEMPT",case,i+1,len(target),f"{n:.9f}",flush=True)
        pdir=dest/f"n_{i:03d}"
        entry=dict(n_fm3=float(n),status="NOT_RUN",attempts=[])
        try:
            anchor=choose_anchor(source,float(n))
            nn,np_,radius=interpolate_initial(anchor,128,float(n))
            parent=dict(kind="H_PREVIOUS_ACCEPTED_PROFILE_N128",
                        source_path=anchor["_source_path"],source_sha256=anchor["_source_sha256"],
                        source_n_fm3=anchor["n_b"],source_branch=spec["branch"])
            first=run_state(EXE,spec["model"],float(n),128,3,"droplet",radius,
                      work_dir=pdir/"primary",initial_nn=nn,initial_np=np_,
                      parent=parent,maxiter=1800,ftol=1e-14,timeout=180)
            attempts=[first]
            if not first.get("accepted_stationary") and "profile" in first:
                warm_nn=np.asarray(first["profile"],dtype=float)[:,0]
                warm_np=np.asarray(first["profile"],dtype=float)[:,1]
                if np.all(np.isfinite(warm_nn+warm_np)) and np.min(warm_nn)>=0 and np.min(warm_np)>=0:
                    second=run_state(EXE,spec["model"],float(n),128,3,"droplet",
                      first["radius"] if first.get("radius",0)>0 else radius,
                      work_dir=pdir/"warm_one",initial_nn=warm_nn,initial_np=warm_np,
                      parent=dict(kind="single_reuse_of_returned_finite_profile",path=first["state_path"],
                                  source_sha256=sha(first["state_path"]),
                                  unchanged_ftol=1e-14,unchanged_physical_acceptance=True),
                      maxiter=1800,ftol=1e-14,timeout=180)
                    attempts.append(second)
            accepted=next((q for q in attempts if q.get("accepted_stationary") and q.get("geometry")=="droplet"),None)
            entry.update(status="ACCEPTED_NATIVE_STATIONARY" if accepted else "REJECTED_NATIVE_STATE",
                accepted_path=str(Path(accepted["state_path"]).relative_to(ROOT)) if accepted else None,
                accepted_sha256=sha(accepted["state_path"]) if accepted else None,
                source_anchor=parent,acceptance_source="spatial TF native stationarity and strict constraints",
                attempts=[dict(state_path=str(Path(a["state_path"]).relative_to(ROOT)),
                               sha256=sha(a["state_path"]),actual_exit=a.get("actual_exit_code"),
                               status=a["status"],accepted_stationary=a["accepted_stationary"],
                               optimizer_mode=a.get("optimizer_mode"),kkt=a.get("kkt"),
                               wall_seconds=a.get("wall_seconds")) for a in attempts])
        except Exception as exc:
            entry.update(status="ENVIRONMENT_OR_INITIALIZATION_FAILURE",error_type=type(exc).__name__,
                         error=str(exc),traceback=traceback.format_exc()[-1700:])
        entry["wall_seconds"]=time.monotonic()-t
        rec["records"].append(entry)
        save(dest/"RESULT.json",rec)
        print("TF_RESULT",case,entry["status"],f"{entry['wall_seconds']:.2f}",flush=True)
    rec.update(status="COMPLETED_WITH_FINITE_REJECTIONS",
               accepted_count=sum(q["status"]=="ACCEPTED_NATIVE_STATIONARY" for q in rec["records"]),
               rejected_count=sum(q["status"]=="REJECTED_NATIVE_STATE" for q in rec["records"]),
               error_count=sum(q["status"]=="ENVIRONMENT_OR_INITIALIZATION_FAILURE" for q in rec["records"]))
    if rec["accepted_count"]==len(target):rec["status"]="ALL_NATIVE_STATIONARY_ACCEPTED"
    save(dest/"RESULT.json",rec)
    print("TF_SUMMARY",case,rec["status"],rec["accepted_count"],len(target),flush=True)
    return rec
if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--case",choices=list(CASES),required=True)
    args=parser.parse_args()
    r=compute(args.case)
    raise SystemExit(0 if r["error_count"]==0 else 1)
