#!/usr/bin/env python3
"""Actual 1.4 Msun phasewise TOV and Cowling for scoped accepted 3-layer EOSs.

No global maximum-mass or radial instability criterion inferred from a single
central-density bracket; record local dM/dn > 0 as a bounded diagnostic.
"""
from __future__ import annotations
import hashlib,json,sys,time,traceback
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import (NativeSharpCowlingDVODE,
                                              solve_sharp_background_native)
from compose_eos.stellar_cowling import CowlingSettings

DATA=ROOT/"data/cowling_maxwell_native_20261119"
OUT=DATA/"fixed_mass_1p4"
TARGET=1.4
CASES={
  "same_fsu2_maxwell":DATA/"three_layer_original_fsu2/maxwell/REPRESENTATION.json",
  "same_fsu2_css":DATA/"three_layer_original_fsu2/css/REPRESENTATION.json",
  "official_fsu2_css":DATA/"official_unified_css_sharp/FSU2/REPRESENTATION.json",
  "official_tm1e_css":DATA/"official_unified_css_sharp/TM1e/REPRESENTATION.json",
  "official_ddme2_css":DATA/"official_unified_css_sharp/DDME2/REPRESENTATION.json",
  "bbp_bps_fsu2_css":DATA/"bbp_bps/whole_models/GPPVA_FSU2/css/REPRESENTATION.json",
  "bbp_bps_tm1e_css":DATA/"bbp_bps/whole_models/GPPVA_TM1e/css/REPRESENTATION.json",
  "bbp_bps_ddme2_css":DATA/"bbp_bps/whole_models/GPPVA_DDME2/css/REPRESENTATION.json",
  "bbp_bps_gm1_css":DATA/"bbp_bps/whole_models/GM_GM1/css/REPRESENTATION.json",
  "bbp_bps_pk1_css":DATA/"bbp_bps/whole_models/PK1_CANONICAL/css/REPRESENTATION.json",
}
TOV=dict(surface_h=3e-8,rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=4097)
COWLING=CowlingSettings(surface_h=3e-8,rtol=3e-9,atol=3e-12,max_step_h=.001)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+"\n")

def run_one(name,path):
    whole=WholeCrustBarotrope.from_contract(json.loads(path.read_text()))
    eos=SharpContractEnthalpyEOS(whole)
    cap=min(.45,float(eos.density_domain[1])-.01)
    low=.25
    if cap<=low:raise ValueError("EOS has inadequate central density support for 1.4 Msun")
    start=time.perf_counter()
    def one(n):return solve_sharp_background_native(eos,float(n),**TOV)
    lo,hi=one(low),one(cap)
    row=dict(name=name,source_path=str(path.relative_to(ROOT)),source_sha256=sha(path),
             central_bracket=[low,cap],bracket_M_Msun=[lo.mass_msun,hi.mass_msun],
             status="RUNNING")
    if not lo.mass_msun<TARGET<hi.mass_msun:
        row["status"]="UNBRACKETED_TARGET_ON_TESTED_CENTRAL_RANGE"
        return row
    nc=brentq(lambda n:one(n).mass_msun-TARGET,low,cap,xtol=2e-8)
    bg=one(nc)
    if abs(bg.mass_msun-TARGET)>3e-6:
        raise RuntimeError("mass target numerical error too large")
    step=.002
    mleft=one(nc-step).mass_msun;mright=one(nc+step).mass_msun
    positive_slope=(mright-mleft)/(2*step)
    row["native_TOV"]=dict(nc_fm3=float(nc),M_Msun=float(bg.mass_msun),
         R_km=float(bg.radius_km),compactness=float(bg.mass_geom_km/bg.radius_km),
         interface_radius_km=[z["r_km"] for z in bg.interfaces],
         source_r_m_continuity=[z["r_jump_km"]==0 and z["m_jump_km"]==0 for z in bg.interfaces],
         local_dM_dnc_MSun_fm3=float(positive_slope),
         local_ascending_M_branch=bool(positive_slope>0))
    if not positive_slope>0:
        row["status"]="TOV_TARGET_FOUND_BUT_LOCAL_ASCENDING_BRANCH_FAILED"
        return row
    solver=NativeSharpCowlingDVODE(bg,COWLING,samples=4097)
    modes={}
    for conv in ("slow","rapid"):
        roots={}
        scan,trace=solver.scan_brackets(np.linspace(900,12000,94),conversion=conv)
        failures=[]
        for br in scan:
            try:
                r=solver.solve_mode(br,conversion=conv,profile_points=180)
                nodes=(r["W_nodes_per_phase"],r["V_nodes_per_phase"])
                if nodes==(0,0):label="f"
                elif nodes==(1,1):label="p1"
                else:continue
                if label in roots:continue
                roots[label]=dict(frequency_hz=float(r["frequency_hz"]),
                     matching_residual=float(r["matching_residual"]),
                     join_relative=float(r["relative_join_mismatch"]),
                     max_interface_traction=float(max(z["traction_relative_residual"] for z in r["interface_junctions"])),
                     nodes=nodes)
                if len(roots)==2:break
            except Exception as exc:
                failures.append(dict(frequency_bracket=br,type=type(exc).__name__,
                                     error=str(exc)[:240]))
        modes[conv]=dict(status="F_P1_COMPLETE" if set(roots)=={"f","p1"}
                         else "PARTIAL_NODE_RESOLVED_SPECTRUM",
                         accepted_modes=roots,failures=failures[:5],
                         brackets_tested=len(scan))
    row["native_cowling"]=modes
    row["max_interface_residual"]=float(solver.max_interface_residual)
    row["status"]=("FIXED_1P4_NATIVE_TOV_F_P1_BOTH_CONVERSIONS_VERIFIED" if
        all(modes[c]["status"]=="F_P1_COMPLETE" for c in ("slow","rapid"))
        else "FIXED_1P4_NATIVE_TOV_PARTIAL_MODE_IDENTITY")
    row["elapsed_seconds"]=time.perf_counter()-start
    return row

def main():
    out=dict(schema="ns_cowling.fixed_mass_1p4_phasewise.v1",
       reference_mass_MSun=TARGET,
       physics="one conserved baryon charge local beta equilibrium, equilibrium barotropic Cowling",
       metric="-+++",
       conditions=dict(TOV=TOV,mode=vars(COWLING)),
       claims_not_made="global maximum mass/radial stability, true microscopic source unification, full metric QNMs",
       cases={},status="RUNNING")
    for name,path in CASES.items():
        try:
            result=run_one(name,path)
        except Exception as e:
            result=dict(status="CASE_FAILED",type=type(e).__name__,error=str(e),
                        traceback=traceback.format_exc()[-1600:])
        out["cases"][name]=result
        print("MASS_1P4",name,result["status"],
              result.get("native_TOV",{}).get("nc_fm3"),
              result.get("native_TOV",{}).get("R_km"),flush=True)
        save(OUT/"PARTIAL.json",out)
    out["status"]="ALL_CASES_ATTEMPTED"
    save(OUT/"RESULT.json",out)

if __name__=="__main__":main()
