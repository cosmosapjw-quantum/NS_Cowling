#!/usr/bin/env python3
"""C17 finite, single-midpoint local density-branch bridge experiment.

Inspect selected already-rejected C17 density edges without changing
original Maxwell node, eigenroot or shape acceptance thresholds.
The first endpoint, midpoint and last endpoint are independently
reconstructed from original EOS. All numerical mode-solving,
frequency scan, shape correlation and uniqueness decisions in Fortran.

Max one inserted midpoint per named edge, no recursive gate inflation.
Original 40-star benchmark results remain immutable.
"""
from __future__ import annotations
import hashlib,io,json,sys,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
sys.path.insert(0,str(ROOT/"scripts"))
import run_c17_extended_density as c17
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS

OUTPUT=ROOT/"data/cowling_c17_extended_branch_20261122/BOUNDED_DENSITY_MIDPOINTS.json"
ORIGINAL=ROOT/"data/cowling_c17_extended_branch_20261122/REAL_P2_INITIAL"
EDGES=(
    ("M2","p1",.60,.55),
    ("M3","p1",.90,.85),
    ("M3","p1",1.02,1.05),
    ("M3","f",.60,.55),
)
def load_eos(zip_bytes,model):
    spec=c17.MODELS[model]
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as full:
      with zipfile.ZipFile(io.BytesIO(full.read("05_SPG5_PK1_10_WHOLE_STARS.zip"))) as inside:
        blob=inside.read(f"FRESH_5x2_WHOLE_STARS/{spec['model']}/{spec['kind']}/REPRESENTATION.json")
    return SharpContractEnthalpyEOS(WholeCrustBarotrope.from_contract(json.loads(blob)))

def run():
    original=c17.SOURCE.read_bytes()
    if hashlib.sha256(original).hexdigest()!=c17.SHA:
        raise RuntimeError("original source lock drift")
    eos={m:load_eos(original,m) for m in ("M2","M3")}
    canonical={m:json.loads((ORIGINAL/f"{m}.json").read_text()) for m in ("M2","M3")}
    rows=[]
    cache={}
    for model,mode,start_factor,end_factor in EDGES:
        label=(model,mode,start_factor,end_factor)
        candidates={}
        for factor in (start_factor,(start_factor+end_factor)/2,end_factor):
            key=(model,round(factor,8))
            if key not in cache:
                try:
                    packet,timing=c17.candidates_for_density(
                        eos[model],c17.MODELS[model]["nc"]*factor,
                        c17.MODELS[model])
                    cache[key]=dict(ok=True,packet=packet)
                except (ValueError,RuntimeError) as e:
                    cache[key]=dict(ok=False,reason=str(e))
            candidates[factor]=cache[key]
        row=dict(model=model,mode=mode,
                 original_status="C17_ORIGINAL_BRANCH_GATE_UNCHANGED",
                 factors=[start_factor,(start_factor+end_factor)/2,end_factor])
        target=0 if mode=="f" else 1
        if not all(v["ok"] for v in candidates.values()):
            row.update(status="BACKGROUND_OR_SCAN_REJECTED",
               errors={str(f):v.get("reason") for f,v in candidates.items() if not v["ok"]})
        else:
            initial=candidates[start_factor]["packet"]["candidate_eigenvalues"]
            original_root=next(r for r in canonical[model]["densities"]
              if r["density_factor"]==start_factor)["local_branches"][mode]["frequency_hz"]
            anchor=c17.select_branch_native(initial,original_root,target)
            if anchor["status"]!="LOCAL_NODE_CONTINUATION_ADMITTED":
                row.update(status="ANCHOR_NOT_RECOVERED",
                    original_anchor_hz=original_root,anchor_result=anchor)
            else:
                first_profile=initial[anchor["candidate_index"]].get("_raw_physical_profile")
                mid_factor=(start_factor+end_factor)/2
                middle=candidates[mid_factor]["packet"]["candidate_eigenvalues"]
                mid=c17.select_branch_with_shape_native(middle,
                       anchor["frequency_hz"],first_profile,target)
                row["first_hop"]=mid
                row["original_anchor_hz"]=original_root
                if mid["status"]!="OVERLAP_AND_NODE_CONTINUATION_ADMITTED":
                    row.update(status="MIDPOINT_GATE_BLOCKED",
                               second_hop=dict(status="NOT_ATTEMPTED_AFTER_MIDPOINT_FAILURE"))
                else:
                    mprofile=middle[mid["candidate_index"]].get("_raw_physical_profile")
                    end=c17.select_branch_with_shape_native(
                       candidates[end_factor]["packet"]["candidate_eigenvalues"],
                       mid["frequency_hz"],mprofile,target)
                    row["second_hop"]=end
                    row["status"]=("TWO_HOP_FORTRAN_CONTINUATION_ADMITTED"
                       if end["status"]=="OVERLAP_AND_NODE_CONTINUATION_ADMITTED"
                       else "ENDPOINT_GATE_STILL_BLOCKED")
        # Drop transient full physical profiles from the persisted JSON
        # and record only bounded outcomes, preserving failure evidence.
        rows.append(row)
        print("C17_MIDPOINT",label,row["status"],
             row.get("first_hop",{}).get("frequency_hz"),
             row.get("second_hop",{}).get("frequency_hz"),flush=True)
    out=dict(schema="ns_cowling.c17.one_midpoint_density_refinement.v1",
             original_source_sha256=c17.SHA,
             original_C17_branch_results="REAL_P2_INITIAL/M2.json + M3.json unchanged",
             source_model_and_node_matching="native Fortran original cs2 + slow/rapid Maxwell + phasewise overlap",
             max_midpoints_per_edge=1,
             physical_acceptance_change=False,
             results=rows)
    OUTPUT.write_text(json.dumps(out,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("C17_MIDPOINT_SAVED",OUTPUT,flush=True)
if __name__=="__main__":run()
