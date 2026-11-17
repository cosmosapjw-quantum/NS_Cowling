#!/usr/bin/env python3
"""H: same-model N128 refinement and actual low-density continuation of two G branches."""
from pathlib import Path
import argparse,hashlib,json,os,sys,time,datetime
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
import repair_maxwell_support_20261117 as g
from compose_eos.tf_candidate import sha,save
BASE=ROOT/"data/low_density_crust_20261117/support_refinement"
GBASE=ROOT/"data/crust_refinement_20261117/support_repair"
EXE=Path("/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117/build/tf_candidate_G_final/tf_candidate_driver")
MATCH=Path("/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117/build/G_verified_general/compose_eos_fortran")
EXE_SHA="21f00ecd9cd492e37329e5690d342842c19025d4b649f97f99e84e592522c189"
MATCH_SHA="533d55ba2b5264e3e5248395baefc9463b4732b818e3e8af14beb774779f5356"
CONFIGS=[
 dict(id="GPPVA_NL3wrL55",model="GPPVA_NL3wrL55",core_id="new_native:PK1_CANONICAL:fine",
      mn=939.,mp=939.,G_anchor="lower_adjacent_branch",anchor_density=.012,
      transition_grid=np.unique(np.r_[np.linspace(.010,.015,21),[.01201,.01202,.01204,.01208]]).tolist(),
      lower_grid=[.0003,.0005,.00075,.001,.0015,.002,.003,.004,.006,.008,.009]),
 dict(id="GPPVA_TM1e",model="GPPVA_TM1e",core_id="new_native:PK1_CANONICAL:fine",
      mn=938.,mp=938.,G_anchor="upper_F_component",anchor_density=.008,
      transition_grid=np.linspace(.006,.012,25).tolist(),
      lower_grid=[.0003,.0005,.00075,.001,.0015,.002,.003,.004,.005]),
]
SOURCE_FILES=[Path(__file__),ROOT/"scripts/repair_maxwell_support_20261117.py",
 ROOT/"scripts/run_new_model_maxwell_20261117.py",
 ROOT/"python_frontend/compose_eos/tf_candidate.py",
 ROOT/"fortran_eos/src/compose_eos_tf_models.f90",
 ROOT/"fortran_eos/src/compose_eos_spatial_tf.f90",
 ROOT/"fortran_eos/tools/tf_candidate_driver.f90"]
def load(p):return json.loads(Path(p).read_text())
def identity():return {str(p.relative_to(ROOT)):sha(p) for p in SOURCE_FILES}
def seed_for(config):
    run=next(r for r in load(GBASE/"REGENERATION.json")["runs"]
       if r["model"]==config["model"] and r["anchor"]==config["G_anchor"] and r["points"]==64)
    point=min([p for p in run["points_results"] if p["accepted"]],key=lambda p:abs(p["n_b"]-config["anchor_density"]))
    old=load(point["state_path"])
    assert old["accepted_stationary"] and old["points"]==64 and old["executable_sha256"]==EXE_SHA
    meta=dict(base_G_commit="905b845ad0f14d2f1e7f67f9eb802691144850cd",
       selected_G_anchor=config["G_anchor"],source_state=point["state_path"],
       source_state_sha256=sha(point["state_path"]),source_density=old["n_b"],
       role="Exact G accepted profile used only as initial guess for all new N128 calculations")
    return old,meta
def prepare():
    BASE.mkdir(parents=True,exist_ok=True)
    (BASE/"seeds").mkdir(parents=True,exist_ok=True)
    assert sha(EXE)==EXE_SHA and sha(MATCH)==MATCH_SHA
    plan=dict(schema="ns_cowling.H_support_refinement_plan.v1",
       base_G_commit="905b845ad0f14d2f1e7f67f9eb802691144850cd",
       configurations=CONFIGS,spatial_mesh=128,source_sha256=identity(),
       TF_executable=str(EXE),TF_executable_sha256=EXE_SHA,matching_executable=str(MATCH),
       matching_executable_sha256=MATCH_SHA,
       units=dict(n="fm^-3",epsilon_pressure="MeV fm^-3",mu="MeV",radius="fm"),
       branch_policy="Continue each selected G branch separately at N128. Initial profile interpolation creates no physical source rows. Failures split support.",
       repair_policy="One eligible optimizer warm restart from its finite profile; unchanged ftol=1e-14, maxiter1800, native physical and optimizer gates.",
       global_minimum_claim=False,continuum_limit_claim=False)
    if (BASE/"PLAN.json").exists():assert load(BASE/"PLAN.json")==plan
    else:save(BASE/"PLAN.json",plan)
    for config in CONFIGS:
        old,meta=seed_for(config);save(BASE/"seeds"/(config["model"]+".json"),meta)
    print(json.dumps(dict(action="prepare",states=sum(len(c["transition_grid"])+len(c["lower_grid"]) for c in CONFIGS),
        plan_sha256=sha(BASE/"PLAN.json"),source_sha256=identity())),flush=True)
    return plan
def generate():
    plan=prepare();runs=[];started=time.monotonic()
    for config in CONFIGS:
        original,seed=seed_for(config)
        grid=np.unique(np.r_[config["lower_grid"],config["transition_grid"]])
        ai=int(np.argmin(abs(grid-config["anchor_density"])))
        anchor="continued_G_"+config["G_anchor"]
        destination=BASE/"states"/config["model"]/anchor/"N128"
        records=[None]*len(grid);last=original;anchor_state=None
        order=[ai]+list(range(ai+1,len(grid)))+list(range(ai-1,-1,-1))
        for index in order:
            if index==ai-1:last=anchor_state or original
            parent=dict(seed=str(BASE/"seeds"/(config["model"]+".json")),
                seed_sha256=sha(BASE/"seeds"/(config["model"]+".json")),
                previous_state=last.get("state_path"),target_mesh=128,
                initial_source="same frozen physical branch; initial guesses are reoptimized natively")
            chosen,point=g.solve_point(EXE,config,grid[index],128,last,seed,destination/f"point_{index:03d}",parent)
            point["grid_index"]=index;records[index]=point
            if chosen is not None:last=chosen
            if index==ai:anchor_state=chosen
            save(BASE/"LIVE_PROGRESS.json",dict(model=config["model"],n_b=float(grid[index]),
                index=index,accepted=point["accepted"],attempts=len(point["attempts"]),
                completed_points=sum(p is not None for p in records),requested_points=len(grid),
                wall_seconds=time.monotonic()-started))
            print(json.dumps(dict(model=config["model"],n_b=float(grid[index]),accepted=point["accepted"],
                                 attempts=len(point["attempts"]))),flush=True)
        record=dict(model=config["model"],config_id=config["id"],anchor=anchor,points=128,
             density_grid_fm3=grid.tolist(),points_results=records,seed=seed,
             accepted_points=sum(p["accepted"] for p in records),requested_points=len(grid))
        save(destination/"RUN.json",record);runs.append(record)
        save(BASE/"REGENERATION_PROGRESS.json",dict(runs=runs))
    assert identity()==plan["source_sha256"]
    result=dict(schema="ns_cowling.H_support_regeneration.v1",plan_sha256=sha(BASE/"PLAN.json"),
       source_sha256=identity(),executable_sha256=EXE_SHA,runs=runs,
       wall_seconds=time.monotonic()-started,
       actual_native_invocations=sum(len(p["attempts"]) for r in runs for p in r["points_results"]),
       accepted_grid_states=sum(r["accepted_points"] for r in runs),
       requested_grid_states=sum(r["requested_points"] for r in runs),source_bytes_unchanged_during_production=True)
    save(BASE/"REGENERATION.json",result)
    print(json.dumps({k:v for k,v in result.items() if k not in ("runs","source_sha256")}),flush=True)
def match():
    os.environ["COMPOSE_EOS_FORTRAN_ROOT"]=str(ROOT)
    os.environ["COMPOSE_EOS_FORTRAN_EXE"]=str(MATCH)
    g.BASE=BASE;g.MATCH_EXE_SHA=MATCH_SHA;g.CONFIGS=CONFIGS
    g.match()
    result=load(BASE/"MATCHING.json")
    result["H_provenance"]=dict(base_G_commit="905b845ad0f14d2f1e7f67f9eb802691144850cd",
         delegated_matching_implementation="scripts/repair_maxwell_support_20261117.py",
         implementation_sha256=sha(ROOT/"scripts/repair_maxwell_support_20261117.py"),
         no_shared_source_mutation=True)
    save(BASE/"MATCHING.json",result)
if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("action",choices=["prepare","generate","match","all"])
    args=parser.parse_args()
    if args.action=="prepare":prepare()
    if args.action in ("generate","all"):generate()
    if args.action in ("match","all"):match()
