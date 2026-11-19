#!/usr/bin/env python3
"""Build fixed 1.4 Msun native Cowling stars with phase-local first-law correction.

Compare to earlier cross-seam pressure-primary closure at the SAME mass,
not just at equal central density. Source original inner/core phase split
is retained. Differences are physical EOS-model choices, not integrator errors.
"""
from __future__ import annotations
import hashlib,json,sys,traceback
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from run_fixed_mass_maxwell_css_20261119 import run_one,TOV,COWLING
BASE=ROOT/"data/cowling_eos_source_calibration_20261119"
M03=ROOT/"data/cowling_maxwell_native_20261119"
OUT=BASE/"phase_local_fixed_mass"
ORIGINAL=json.loads((M03/"fixed_mass_1p4/RESULT.json").read_text())
SOURCE=json.loads((BASE/"phase_local_repairs/RESULT.json").read_text())

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+"\n")

def one(name):
    p=BASE/"phase_local_repairs/whole"/name/"phasewise_firstlaw_integrated/css/REPRESENTATION.json"
    result=run_one(name,p)
    result["source_phase_local_repair"]=True
    result["source_contract_sha256"]=sha(p)
    result["source_case_hash"]=sha(BASE/"phase_local_repairs/RESULT.json")
    if result["status"]=="FIXED_1P4_NATIVE_TOV_F_P1_BOTH_CONVERSIONS_VERIFIED":
        older=ORIGINAL["cases"]["official_"+name.lower()+"_css"]
        if older["status"]!="FIXED_1P4_NATIVE_TOV_F_P1_BOTH_CONVERSIONS_VERIFIED":
            raise RuntimeError("previous original unified CSS control is not source-locked")
        curr=result["native_TOV"];prev=older["native_TOV"]
        diff=dict(M_msun=curr["M_Msun"]-prev["M_Msun"],
                  R_km=curr["R_km"]-prev["R_km"],
                  central_density_fm3=curr["nc_fm3"]-prev["nc_fm3"],
                  mode_frequency_hz={})
        for conv in ("slow","rapid"):
            diff["mode_frequency_hz"][conv]={
                label:result["native_cowling"][conv]["accepted_modes"][label]["frequency_hz"]-
                      older["native_cowling"][conv]["accepted_modes"][label]["frequency_hz"]
                for label in ("f","p1")}
        result["phasewise_minus_cross_seam_firstlaw"]=diff
        result["comparison_physics_warning"]="Different Euler-consistent phase bridge; NOT a numerical solver error or original unaltered source EOS"
    return result

def main():
    x=dict(schema="ns_cowling.phase_local_fixed_1p4_comparison.v1",
      reference_old_cross_seam_result_sha256=sha(M03/"fixed_mass_1p4/RESULT.json"),
      source_repair_full_evidence_sha256=sha(BASE/"phase_local_repairs/RESULT.json"),
      same_mass_MSun=1.4,
      methods="native Fortran phasewise TOV and matched slow/rapid Cowling",
      source_claim="original P(n) source knots unchanged, epsilon(mu) thermodynamically repaired independently within original Iphase7 and Iphase6, original phase-gap never integrated over",
      cases={},status="RUNNING")
    for name in ("FSU2","TM1e","DDME2"):
        try:row=one(name)
        except Exception as ex:row=dict(status="CASE_FAILED",type=type(ex).__name__,
             error=str(ex),traceback=traceback.format_exc()[-1400:])
        x["cases"][name]=row
        print("PHASE_LOCAL_FIXED_MASS",name,row["status"],
              row.get("native_TOV",{}).get("R_km"),
              row.get("phasewise_minus_cross_seam_firstlaw",{}).get("R_km"),
              flush=True)
        save(OUT/"PARTIAL.json",x)
    x["status"]="ALL_CALIBRATION_CASES_ATTEMPTED"
    save(OUT/"RESULT.json",x)
if __name__=="__main__":main()
