#!/usr/bin/env python3
"""Track the W1 overtone through a node-identity change in phase-local TM1e.

This numerical experiment records individual W and V nodes. A high shape
overlap is evidence of continuity, but a changed Sturm-node pair is not
silently relabelled p1 without a physical mode-family theorem.
"""
from __future__ import annotations
import hashlib,json,sys,time,traceback
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import (NativeSharpCowlingDVODE,solve_sharp_background_native)
from compose_eos.stellar_mode_continuation import eigenfunction_overlap,track_mode_family
from compose_eos.stellar_cowling import CowlingSettings

BASE=ROOT/"data/cowling_eos_source_calibration_20261119"
SOURCE=BASE/"phase_local_repairs/whole/TM1e/phasewise_firstlaw_integrated/css/REPRESENTATION.json"
OUT=BASE/"tm1e_mode_identity"
DENSITIES=(.30,.32,.34,.35,.36,.36680854236573207,.38)
TOV=dict(surface_h=3e-8,rtol=3e-10,atol=3e-11,
         max_step_h=.004,native_samples=8193)
MODES=CowlingSettings(surface_h=3e-8,rtol=3e-10,atol=3e-13,max_step_h=.0005)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,r):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(r,indent=2,ensure_ascii=False,allow_nan=False)+"\n")

def run():
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(SOURCE.read_text())))
    out=dict(schema="ns_cowling.tm1e_phase_local_W1_overtone_node_change.v1",
       source_sha256=sha(SOURCE),source_path=str(SOURCE.relative_to(ROOT)),
       density_list=list(DENSITIES),params=dict(TOV=TOV,Cowling=vars(MODES)),
       physical_claim="candidate W1 overtone, label p1 withheld if V1 node condition fails",
       status="RUNNING",cases={})
    for conv in ("slow","rapid"):
        previous=None;trace=[]
        for nc in DENSITIES:
            bg=solve_sharp_background_native(eos,nc,**TOV)
            sol=NativeSharpCowlingDVODE(bg,MODES,samples=8193)
            brs,_=sol.scan_brackets(np.linspace(1400,8800,112),conversion=conv)
            chosen=None;accepted=[];reject=[]
            for bracket in brs:
                try:
                    z=sol.solve_mode(bracket,conversion=conv,profile_points=450)
                    accepted.append((float(z["frequency_hz"]),
                                     z["W_nodes_per_phase"],z["V_nodes_per_phase"]))
                    if z["W_nodes_per_phase"]==1 and chosen is None:
                        chosen=z
                except Exception as ex:
                    reject.append(dict(bracket=list(bracket),failure=str(ex)[:200]))
                if chosen is not None:
                    break
            continuity=None
            if chosen is not None and previous is not None:
                shape=eigenfunction_overlap(previous,chosen)
                tracked=track_mode_family([previous],[chosen],min_overlap=.88)
                conditional=track_mode_family([previous],[chosen],min_overlap=.88,
                                        v_node_policy="shape_continuity")
                continuity=dict(physical_profile_shape_overlap=shape,
                                node_pair_safe_tracker=tracked,
                                conditional_shape_continuity_tracker=conditional,
                                previous_nodes=[int(previous["W_nodes_per_phase"]),
                                                int(previous["V_nodes_per_phase"])])
            record=dict(nc_fm3=float(nc),M_MSun=bg.mass_msun,R_km=bg.radius_km,
                status="W1_ROOT_FOUND" if chosen is not None else "NO_W1_MODE_FOUND",
                accepted_first_roots=accepted[:5],
                failed_root_checks=reject[:3],continuation=continuity)
            if chosen is not None:
                record["candidate"]=dict(frequency_hz=float(chosen["frequency_hz"]),
                    W_nodes=int(chosen["W_nodes_per_phase"]),
                    V_nodes=int(chosen["V_nodes_per_phase"]),
                    join_relative=float(chosen["relative_join_mismatch"]),
                    interface_traction_max=float(max(x["traction_relative_residual"]
                                                    for x in chosen["interface_junctions"])))
            previous=chosen
            trace.append(record)
            print("TM1e_OVERTONE",conv,nc,record.get("candidate"),
                  "continuation",continuity,flush=True)
        out["cases"][conv]=trace
        save(OUT/"PARTIAL.json",out)
    out["status"]="NUMERICALLY_TRACED_WITH_NODE_AMBIGUITY_RETAINED"
    save(OUT/"RESULT.json",out)
if __name__=="__main__":run()
