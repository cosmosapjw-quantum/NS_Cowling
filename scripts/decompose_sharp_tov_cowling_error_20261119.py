#!/usr/bin/env python3
"""Separate TOV background error and Cowling coefficient sampling error.

At identical EOS and source-locked brackets:
 A) DOP853 Python TOV + native Cowling, increasingly fine samples;
 B) Native TOV at different samples + fixed 32769 Cowling samples;
 C) Native TOV and native Cowling at matched samples.
No claim of formal asymptotic convergence when oscillatory errors persist.
"""
from __future__ import annotations
import hashlib,json,sys,time,traceback
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS,solve_sharp_background
from compose_eos.stellar_sharp_native import (solve_sharp_background_native,NativeSharpCowlingDVODE)
from compose_eos.stellar_cowling import CowlingSettings

BASE=ROOT/"data/cowling_maxwell_native_20261119"
OUT=ROOT/"data/cowling_eos_source_calibration_20261119/resolution_separation"
PREVIOUS=json.loads((BASE/"validated_convergence/RESULT.json").read_text())
CASES={
 "same_FSU2_maxwell":BASE/"three_layer_original_fsu2/maxwell/REPRESENTATION.json",
 "original_FSU2_css":BASE/"official_unified_css_sharp/FSU2/REPRESENTATION.json",
}
RESOLUTIONS=(4097,8193,16385,32769)
TOV=dict(surface_h=3e-8,rtol=3e-10,atol=3e-11,max_step_h=.004)
SETTING=CowlingSettings(surface_h=3e-8,rtol=3e-9,atol=3e-12,max_step_h=.001)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+"\n")

def modes(background,reference,ns):
    start=time.monotonic()
    solver=NativeSharpCowlingDVODE(background,SETTING,samples=ns)
    setup=time.monotonic()-start
    modes_out={}
    for conversion in ("slow","rapid"):
        modes_out[conversion]={}
        for label in ("f","p1"):
            f=float(reference[conversion][label]);br=(f*.993,f*1.007)
            start=time.monotonic()
            try:
                result=solver.solve_mode(br,conversion=conversion,profile_points=200)
                if (result["W_nodes_per_phase"],result["V_nodes_per_phase"])!=(
                    (0,0) if label=="f" else (1,1)):
                    raise RuntimeError("nodal mode identity changed")
                modes_out[conversion][label]=dict(status="CHECKED",
                    f_Hz=float(result["frequency_hz"]),
                    joint_residual=float(result["relative_join_mismatch"]),
                    wall_seconds=time.monotonic()-start)
            except Exception as e:
                modes_out[conversion][label]=dict(status="REJECTED",
                    error=str(e)[:300],wall_seconds=time.monotonic()-start)
    return dict(coefficient_sample_count=int(ns),
                coefficient_setup_wall_seconds=setup,
                modes=modes_out)

def run_one(name,path):
    ref=PREVIOUS["cases"][name]["input_reference_frequencies_hz"]
    eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(path.read_text())))
    r=dict(name=name,source_sha256=sha(path),controls=dict(tov=TOV,cowling=vars(SETTING)),
           original_20261119_reference=ref,python_background={},native_background={})
    start=time.monotonic()
    python_bg=solve_sharp_background(eos,.3,**TOV)
    r["python_background"]=dict(M_Msun=python_bg.mass_msun,
        R_km=python_bg.radius_km,seconds=time.monotonic()-start,coefficient_sampling={})
    for ns in RESOLUTIONS:
        r["python_background"]["coefficient_sampling"][str(ns)]=modes(python_bg,ref,ns)
        print("PYTHON_TOV_NATIVE_COWLING",name,ns,
              [(cv,k,round(v.get("f_Hz",float("nan")),5))
               for cv,z in r["python_background"]["coefficient_sampling"][str(ns)]["modes"].items()
               for k,v in z.items()],flush=True)
    for ns in RESOLUTIONS:
        start=time.monotonic()
        bg=solve_sharp_background_native(eos,.3,native_samples=ns,**TOV)
        row=dict(M_Msun=bg.mass_msun,R_km=bg.radius_km,
                 M_minus_python_Msun=bg.mass_msun-python_bg.mass_msun,
                 R_minus_python_km=bg.radius_km-python_bg.radius_km,
                 wall_seconds=time.monotonic()-start)
        row["matched_samples"]=modes(bg,ref,ns)
        row["fine_cowling_samples_32769"]=modes(bg,ref,32769)
        r["native_background"][str(ns)]=row
        print("NATIVE_TOV_COEFF_DECOMPOSE",name,ns,"Mdiff",row["M_minus_python_Msun"],
           "Rdiff",row["R_minus_python_km"],flush=True)
        dump(OUT/"PARTIAL.json",dict(case=r,last_ns=ns))
    return r

def main():
    result=dict(schema="ns_cowling.sharp_error_budget_decomposition.v1",
        controls=dict(tov=TOV,mode=vars(SETTING)),cases={},
        interpretation="vary TOV source and coefficient interpolation separately; always same EOS and mode bracket",
        status="RUNNING")
    for name,path in CASES.items():
        try:
            result["cases"][name]=run_one(name,path)
        except Exception as ex:
            result["cases"][name]=dict(status="BLOCKED",type=type(ex).__name__,
                error=str(ex),traceback=traceback.format_exc()[-1700:])
        dump(OUT/"PARTIAL.json",result)
    result["status"]="COMPLETED_WITH_SCOPED_ERROR_COMPONENTS"
    dump(OUT/"RESULT.json",result)
    print("ERROR_BUDGET_FINISH",[(k,"native_background" in v) for k,v in result["cases"].items()],flush=True)
if __name__=="__main__":main()
