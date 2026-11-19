#!/usr/bin/env python3
"""Follow phasewise native f/p1 along connected cold stellar sequences.

This does not prove global mode tracking at avoided crossings. At each fixed
EOS and conversion, compare adjacent eigenfunctions by a phase-local
sign-independent W/rV shape overlap; retain genuine ambiguity.
"""
from __future__ import annotations
import hashlib,json,sys,time,traceback
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_sharp_native import NativeSharpCowlingDVODE,solve_sharp_background_native
from compose_eos.stellar_mode_continuation import track_mode_family,eigenfunction_overlap
from compose_eos.stellar_cowling import CowlingSettings

BASE=ROOT/"data/cowling_maxwell_native_20261119"
OUT=ROOT/"data/cowling_eos_source_calibration_20261119/mode_continuation"
CASES={
    "same_fsu2_maxwell":BASE/"three_layer_original_fsu2/maxwell/REPRESENTATION.json",
    "same_fsu2_css":BASE/"three_layer_original_fsu2/css/REPRESENTATION.json",
    "official_fsu2_corrected_css":BASE/"official_unified_css_sharp/FSU2/REPRESENTATION.json",
}
DENSITIES=(.28,.30,.33,.36,.39)
TOV=dict(surface_h=3e-8,rtol=3e-10,atol=3e-11,max_step_h=.004,native_samples=4097)
SETTING=CowlingSettings(surface_h=3e-8,rtol=3e-9,atol=3e-12,max_step_h=.001)

def sha(x):return hashlib.sha256(x.read_bytes()).hexdigest()
def dump(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+"\n")

def mode_summary(m):
    return dict(frequency_hz=float(m["frequency_hz"]),W_nodes=int(m["W_nodes_per_phase"]),
                V_nodes=int(m["V_nodes_per_phase"]),relative_join=float(m["relative_join_mismatch"]),
                max_interface_traction_relative=float(max(z["traction_relative_residual"]
                                                for z in m["interface_junctions"])),
                interface_jump_count=len(m["interface_junctions"]),
                root_iterations=int(m["root_iterations"]))

def modes_at_density(eos,n,conversion):
    start=time.monotonic()
    bg=solve_sharp_background_native(eos,n,**TOV)
    solver=NativeSharpCowlingDVODE(bg,SETTING,samples=4097)
    brackets,scan=solver.scan_brackets(np.linspace(1250.,8000.,82),conversion=conversion)
    roots={}
    failures=[]
    for br in brackets:
        try:
            mode=solver.solve_mode(br,conversion=conversion,profile_points=360)
            ident=(mode["W_nodes_per_phase"],mode["V_nodes_per_phase"])
            if ident==(0,0):label="f"
            elif ident==(1,1):label="p1"
            else:continue
            if label not in roots:roots[label]=mode
            if set(roots)=={"f","p1"}:break
        except Exception as ex:
            failures.append(dict(bracket=br,type=type(ex).__name__,
                                 message=str(ex)[:180]))
    return bg,roots,failures,len(brackets),time.monotonic()-start

def main():
    out=dict(schema="ns_cowling.cowling_source_mode_family_continuation.v1",
        density_grid_fm3=list(DENSITIES),controls=dict(TOV=TOV,Cowling=vars(SETTING)),
        overlap_diagnostic="phase-separated equal weighted W/rV shape on normalized local radius",
        no_global_sturm_liouville_theorem_or_avoided_crossing_claim=True,
        cases={},status="RUNNING")
    for name,path in CASES.items():
        entry=dict(source_contract_sha256=sha(path),source_path=str(path.relative_to(ROOT)),
                   conversion={},status="RUNNING")
        out["cases"][name]=entry
        try:
            eos=SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(path.read_text())))
            for conversion in ("slow","rapid"):
                rows=[];previous_modes=None
                for density in DENSITIES:
                    bg,modes,failures,nbrackets,seconds=modes_at_density(eos,density,conversion)
                    labels=[k for k in ("f","p1") if k in modes]
                    selected=[modes[k] for k in labels]
                    identity=None
                    if previous_modes is not None and len(previous_modes)==2 and len(selected)==2:
                        identity=track_mode_family(previous_modes,selected,
                            min_overlap=.88,max_fractional_step=.3)
                    row=dict(n_c_fm3=float(density),M_MSun=float(bg.mass_msun),
                        R_km=float(bg.radius_km),interface_r_km=[x["r_km"] for x in bg.interfaces],
                        frequencies={k:mode_summary(v) for k,v in modes.items()},
                        candidate_brackets=int(nbrackets),seconds=float(seconds),
                        failures=failures[:4],continuation_from_previous=identity,
                        status="F_P1_NODE_IDENTIFIED" if len(modes)==2 else "PARTIAL_SPECTRUM")
                    rows.append(row)
                    print("MODE_FAMILY",name,conversion,density,row["status"],
                          [(k,round(v["frequency_hz"],5)) for k,v in row["frequencies"].items()],
                          "overlap",identity,flush=True)
                    previous_modes=selected if len(selected)==2 else None
                entry["conversion"][conversion]=rows
            good=all(row["status"]=="F_P1_NODE_IDENTIFIED"
                     and (row["continuation_from_previous"] is None
                          or row["continuation_from_previous"]["status"]=="CONTINUED")
                     for rows in entry["conversion"].values() for row in rows)
            entry["status"]="SCOPED_MODE_FAMILY_CONTINUED" if good else "MODE_FAMILY_PARTIALLY_UNRESOLVED"
        except Exception as ex:
            entry["status"]="CASE_BLOCKED";entry["reason"]=dict(
                type=type(ex).__name__,message=str(ex),
                traceback=traceback.format_exc()[-1200:])
        dump(OUT/"PARTIAL.json",out)
    out["status"]="BOUNDED_CONTINUATION_COMPLETE"
    dump(OUT/"RESULT.json",out)
    print("MODE_FAMILY_END",[(k,v["status"]) for k,v in out["cases"].items()],flush=True)

if __name__=="__main__":main()
