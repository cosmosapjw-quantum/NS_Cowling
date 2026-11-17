#!/usr/bin/env python3
"""One diagnosed opposite-neighbor N128 control at a support-blocking state."""
from pathlib import Path
import os,sys,json,copy,time
ROOT=Path("/home/cosmosapjw/ns_cowling_wt/low_density_crust_20261117")
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
import refine_supported_crust_20261117 as h
import repair_maxwell_support_20261117 as g
from compose_eos.tf_candidate import sha,save,interpolate_initial,run_state
BASE=Path(__file__).resolve().parent
def load(p):return json.loads(Path(p).read_text())
def main():
    regeneration=load(h.BASE/"REGENERATION.json")
    run=next(r for r in regeneration["runs"] if r["model"]=="GPPVA_TM1e")
    point=next(p for p in run["points_results"] if abs(p["n_b"]-.008)<1e-12)
    neighbor=next(p for p in run["points_results"] if abs(p["n_b"]-.00825)<1e-12)
    assert not point["accepted"] and neighbor["accepted"]
    source=load(neighbor["state_path"])
    plan=dict(schema="ns_cowling.H_bounded_anchor_control.v1",density_fm3=.008,mesh=128,
       reason="The sole source gap between accepted n=.00775 and .00825 prevents a Maxwell crossing in the connected low-density component.",
       mechanism="Previous seed was interpolated from G N64 at n=.008. This control interpolates the accepted same-branch N128 neighbor n=.00825 downward, preserving all solver gates.",
       source_state=neighbor["state_path"],source_state_sha256=sha(neighbor["state_path"]),
       replaced_failed_point=point,original_regeneration_sha256=sha(h.BASE/"REGENERATION.json"),
       code_sha256=sha(__file__),frozen_API_sha256=sha(ROOT/"python_frontend/compose_eos/tf_candidate.py"),
       executable_sha256=sha(h.EXE),maxiter=1800,ftol=1e-14,max_new_invocations=1,
       no_additional_warm_restart=True,other_failed_points_preserved=[.0085,.0095])
    save(BASE/"PLAN.json",plan)
    started=time.monotonic()
    nn,pp,radius=interpolate_initial(source,128,.008)
    s=run_state(h.EXE,"GPPVA_TM1e",.008,128,3,"droplet",radius,work_dir=BASE/"state",
       initial_nn=nn,initial_np=pp,parent=dict(kind="bounded_opposite_neighbor",source=neighbor["state_path"],
         source_sha256=sha(neighbor["state_path"]),plan=str(BASE/"PLAN.json"),plan_sha256=sha(BASE/"PLAN.json")),
       maxiter=1800,ftol=1e-14)
    old=load(point["attempts"][-1]["state_path"])
    comparison=dict(E_per_baryon_change_MeV=s["epsilon"]/.008-old["epsilon"]/.008,
       radius_change_fm=s["radius"]-old["radius"],
       reason="Comparison diagnoses the finite failed iterate; it does not admit or export that iterate.")
    receipt=dict(accepted=s["accepted_stationary"],state_path=s["state_path"],state_sha256=sha(s["state_path"]),
       actual_native_exit=s["actual_exit_code"],wall_seconds=time.monotonic()-started,comparison=comparison)
    save(BASE/"CONTROL_RESULT.json",receipt)
    if not s["accepted_stationary"]:
        print(json.dumps(dict(status="CONTROL_REJECTED",**receipt)),flush=True);return
    assert abs(comparison["E_per_baryon_change_MeV"])<1e-7
    revised=copy.deepcopy(run)
    replacement=dict(n_b=.008,points=128,model="GPPVA_TM1e",accepted=True,state_path=s["state_path"],
       attempts=[dict(state_path=s["state_path"],state_sha256=sha(s["state_path"]),status=s["status"],
          accepted_stationary=s["accepted_stationary"],actual_exit_code=s["actual_exit_code"])],
       grid_index=point["grid_index"],supersedes_failed_point_evidence=point,control_plan_sha256=sha(BASE/"PLAN.json"))
    revised["points_results"][point["grid_index"]]=replacement
    revised["accepted_points"]=sum(p["accepted"] for p in revised["points_results"])
    revised["original_run_sha256"]=sha(h.BASE/"states/GPPVA_TM1e/continued_G_upper_F_component/N128/RUN.json")
    out=dict(schema="ns_cowling.H_bounded_control_regeneration.v1",runs=[revised],
       original_regeneration_sha256=sha(h.BASE/"REGENERATION.json"),control=receipt,actual_new_native_invocations=1)
    save(BASE/"REGENERATION.json",out)
    os.environ["COMPOSE_EOS_FORTRAN_ROOT"]=str(ROOT);os.environ["COMPOSE_EOS_FORTRAN_EXE"]=str(h.MATCH)
    g.BASE=BASE;g.MATCH_EXE_SHA=h.MATCH_SHA;g.CONFIGS=h.CONFIGS
    g.match()
    print(json.dumps(dict(status="CONTROL_ACCEPTED",**receipt)),flush=True)
if __name__=="__main__":main()
