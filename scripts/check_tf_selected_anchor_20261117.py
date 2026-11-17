#!/usr/bin/env python3
"""Exactly four DDME2 anchors from the actual F-selected retained profile.

The production droplet-seeded branch remains immutable. Each control independently
interpolates the same historical selected profile, then solves the native fields.
"""
from pathlib import Path
import csv,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend"),str(ROOT/"pipeline")]
from run_crust_refinement_20261117 import solve, BASE, load, sha, save
from compose_eos.tf_candidate import weights
folder=ROOT/"data/research_native_20261116/source_grid_cohort/GPPVA_DDME2/attempt_001/output/generated"
selected=min(csv.DictReader((folder/"selected.csv").open()),key=lambda s:abs(float(s["n_b"])-.01))
attempt=next(s for s in csv.DictReader((folder/"attempts.csv").open())
    if s["node"]==selected["node"] and s["mesh"]==selected["mesh"] and s["attempt"]==selected["attempt"])
raw=[s for s in csv.DictReader((folder/"profiles.csv").open())
    if s["node"]==selected["node"] and s["mesh"]==selected["mesh"] and s["attempt"]==selected["attempt"]]
raw.sort(key=lambda s:int(s["cell"]))
seed={"id":f"legacy_DDME2_selected_node{selected['node']}_attempt{selected['attempt']}",
      "n_b":float(selected["n_b"]),"points":int(selected["mesh"]),"dimension":int(attempt["dimension"]),
      "seed":attempt["seed"],"radius":float(attempt["radius"]),
      "profile":[[float(s[k]) for k in ["nn","np","ne","mu_n","mu_p","mu_e","sigma","omega","rho","V"]] for s in raw],
      "selected_row":selected,"selected_attempt":attempt,
      "source_sha256":{str((folder/n).relative_to(ROOT)):sha(folder/n) for n in ["selected.csv","attempts.csv","profiles.csv"]}}
destination=BASE/"selected_anchor_control"
destination.mkdir(exist_ok=True)
save(destination/"SEED.json",seed)
results=[]
for points in [16,32,64,128]:
    s=solve(seed,"DDME2",.01,points,f"SELECTED_F_ANCHOR_N{points:03d}")
    original=load(BASE/f"DDME2_N{points:03d}_STATES.json")[4]
    p=np.array(s["profile"]); q=np.array(original["profile"]); w=weights(points,3)
    rec={"points":points,"n_b":.01,"control_state":s["state_path"],"control_sha256":sha(s["state_path"]),
         "production_state":original["state_path"],"production_sha256":sha(original["state_path"]),
         "accepted":s["accepted_stationary"],"control_geometry":s["geometry"],"production_geometry":original["geometry"],
         "control_scalars":{k:s[k] for k in ["epsilon","pressure","mu_B","radius","kkt","radius_stationarity"]},
         "delta_control_minus_production":{k:s[k]-original[k] for k in ["epsilon","pressure","mu_B","radius"]},
         "delta_energy_per_baryon":(s["epsilon"]-original["epsilon"])/.01,
         "weighted_rms_nucleon_profile_difference_fm3":float(np.sqrt(w@((p[:,0]-q[:,0])**2+(p[:,1]-q[:,1])**2))),
         "max_nucleon_profile_difference_fm3":float(np.max(np.abs(p[:,:2]-q[:,:2]))),
         "comparison_scope":"same physical model, density, mesh, native executable and unchanged acceptance gates; initial basin differs"}
    results.append(rec)
save(destination/"SUMMARY.json",{"schema":"ns_cowling.tf_selected_anchor_control.v1",
     "script_sha256":sha(Path(__file__)),"seed_sha256":sha(destination/"SEED.json"),"results":results,
     "scope":"four fixed-density controls only; no replacement of production state or full-grid branch-minimum claim"})
print("COMPLETED selected_anchor_control",flush=True)
