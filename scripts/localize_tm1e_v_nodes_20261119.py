#!/usr/bin/env python3
"""Inspect slow TM1e W1 branch V-node turning within each physical phase.

Repeat eigenfunction reconstruction with distinct output grids. This checks
whether the V-node transition in Cowling is a plotting/discretization artifact.
It is not a general mode classification theorem.
"""
from __future__ import annotations

import json,sys,hashlib
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import NativeSharpCowlingDVODE,solve_sharp_background_native
from compose_eos.stellar_cowling import CowlingSettings

SOURCE=ROOT/"data/cowling_eos_source_calibration_20261119/phase_local_repairs/whole/TM1e/phasewise_firstlaw_integrated/css/REPRESENTATION.json"
OUT=ROOT/"data/cowling_eos_source_calibration_20261119/tm1e_mode_identity/NODE_LOCALIZATION.json"
DATA={.34:6553.1706,.35:6623.3696,.36680854236573207:6692.2003,.38:6719.2902}
TOV=dict(surface_h=3e-8,rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=8193)
SET=CowlingSettings(surface_h=3e-8,rtol=3e-10,atol=3e-13,max_step_h=.0005)


def crossings(points,key):
    a=np.array([[p["r_km"],p[key]] for p in points],float)
    r,y=a.T
    scale=max(abs(y))
    if scale==0:return []
    mask=abs(y)>scale*1e-9
    r,y=r[mask],y[mask]
    pts=[]
    for i in np.flatnonzero(y[1:]*y[:-1]<0):
        root=float(r[i]-y[i]*(r[i+1]-r[i])/(y[i+1]-y[i]))
        pts.append(root)
    return pts

def run():
    eo=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(SOURCE.read_text())))
    record=dict(schema="ns_cowling.tm1e_v_node_localization.v1",
      source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
      controls=dict(TOV=TOV,Cowling=vars(SET)),
      phase_contiguous_node_count=True,cases={})
    for nc,expected in DATA.items():
        bg=solve_sharp_background_native(eo,nc,**TOV)
        solver=NativeSharpCowlingDVODE(bg,SET,samples=8193)
        case=dict(central_density_fm3=nc,physical_interface_r_km=bg.interfaces[0]["r_km"],
                  surface_r_km=bg.radius_km,results={})
        for output_n in (160,450,900,1800):
            sol=solver.solve_mode((expected-15,expected+15),
                                 conversion="slow",profile_points=output_n)
            phase_data={}
            for phase in (0,1):
                records=[p for p in sol["profile"] if p["phase_index"]==phase]
                phase_data[str(phase)]=dict(
                    W_crossings_km=crossings(records,"W"),
                    V_crossings_km=crossings(records,"V_km_inverse"),
                    radius_min_km=min(z["r_km"] for z in records),
                    radius_max_km=max(z["r_km"] for z in records),
                    source_profile_points=len(records))
            case["results"][str(output_n)]=dict(
                frequency_hz=sol["frequency_hz"],
                nodes_W=sol["W_nodes_per_phase"],nodes_V=sol["V_nodes_per_phase"],
                phases=phase_data,relative_join=sol["relative_join_mismatch"])
            print("TM1e_NODE_LOCALIZATION",nc,output_n,
                "freq",sol["frequency_hz"],"nodes",(sol["W_nodes_per_phase"],sol["V_nodes_per_phase"]),
                "Wzeros",[(k,v["W_crossings_km"]) for k,v in phase_data.items()],
                "Vzeros",[(k,v["V_crossings_km"]) for k,v in phase_data.items()],flush=True)
        record["cases"][str(nc)]=case
    OUT.write_text(json.dumps(record,indent=2)+"\n")
    print("NODE_LOCALIZATION_RECORDED",OUT,flush=True)

if __name__=="__main__":run()
