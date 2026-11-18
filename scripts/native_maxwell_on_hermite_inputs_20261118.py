"""Run native Maxwell on identical source-tier EOS after cubic/quintic
source interpolation. CSS was run separately on these same representations.

Inputs are the 7 previously verified finite candidates, not all 635
(the others were rejected/blocked in source-level screening). No original
source row is overwritten. Source physical components are restricted to the
documented locally crossing phase branches, never glued across gaps.
"""
from __future__ import annotations
import hashlib,json,time,sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from all_origin_maxwell_css_20261118 import OUT,sha
from compare_hermite_as_interpolation_only_20261118 import (
    resample_eos,source,CAT,INTERPOLATION_METHODS,MATCHING_METHODS,save)
from check_hermite_roots_for_dual_pairs_20261118 import comps
from compose_eos.api import construct_maxwell_transition,write_compose,read_compose
from compose_eos.maxwell_v2 import _segment_eosdata

BASE=OUT/"HERMITED_INPUT_2X2"
DEST=BASE/"NATIVE_MAXWELL"

def component_in_source(data,entry,label):
    c=next((x for x in comps(data,entry,False) if x[2]==label),None)
    if c is None:raise ValueError("specified locally crossing source component not found: "+label)
    return c

def slice_selected(source_data,source_entry,derived,label):
    c=component_in_source(source_data,source_entry,label)
    nlow=float(source_data.nb[c[0]]);nhigh=float(source_data.nb[c[1]-1])
    positions=np.flatnonzero((derived.nb>=nlow-1e-13)&(derived.nb<=nhigh+1e-13))
    if len(positions)<3 or positions[-1]-positions[0]+1!=len(positions):
        raise ValueError("derived interpolant breaks admitted source component contiguity")
    arr=_segment_eosdata(derived,(0,int(positions[0]),int(positions[-1])+1))
    if not (np.isclose(arr.nb[0],nlow,rtol=0,atol=1e-13) and
            np.isclose(arr.nb[-1],nhigh,rtol=0,atol=1e-13)):
        raise ValueError("source component endpoints changed by interpolation")
    return arr,dict(component=label,original_n_min=nlow,original_n_max=nhigh,
                    original_rows=c[1]-c[0],interpolated_rows=arr.n)

def one(pair,order):
    tier,crust,core=pair["tier"],pair["crust"],pair["core"]
    dest=DEST/tier/(crust+"__"+core)/order
    if dest.exists():raise FileExistsError(str(dest))
    dest.mkdir(parents=True)
    case=json.loads((BASE/tier/(crust+"__"+core)/"RESULT.json").read_text())
    root=case["interpolation_pair_results"][order+"+maxwell"]["roots"]
    if len(root)!=1:raise ValueError("require exactly one sealed cubic/quintic local Maxwell crossing")
    result=dict(status="RUNNING",tier=tier,crust=crust,core=core,
        order=order,matching="maxwell",
        source_ids=case["original_source_ids"],
        source_interpolants_are_derived=True,
        never_cross_tier=True,
        raw_source_rows_not_overwritten=True,
        source_component_names=root[0]["source_components"])
    save(dest/"RESULT.json",result)
    try:
        l=source(pair["crust_source_id"]);h=source(pair["core_source_id"])
        dl,rl=resample_eos(l,order)
        dh,rh=resample_eos(h,order)
        low,lc=slice_selected(l,CAT[pair["crust_source_id"]],dl,root[0]["source_components"][0])
        high,hc=slice_selected(h,CAT[pair["core_source_id"]],dh,root[0]["source_components"][1])
        matched,report=construct_maxwell_transition(low,high,
             assume_local_charge_neutrality=True,return_report=True)
        target=dest/"compose"
        write_compose(target,matched)
        back=read_compose(target)
        maxdiff=max(float(np.max(abs(getattr(back,z)-getattr(matched,z))))
                    for z in ("nb","Q1","Q3","Q6","Q7"))
        if maxdiff>1e-9:raise ValueError("native source matched serialization failed")
        p=float(report["pressure_transition"]);mu=float(report["mu_transition"])
        nlo=float(report["n_low"]);nhi=float(report["n_high"])
        latent=float(report["latent_energy_density"])
        if not(0<nlo<nhi and p>0 and mu>0 and latent>0):
            raise ValueError("native Maxwell physical support gate")
        identity=abs(latent-mu*(nhi-nlo))
        if identity>1e-8:raise ValueError("native Maxwell identity gate")
        result.update(status="PASS_NATIVE_MAXWELL_MATCHED",
            transition_P=p,mu_t=mu,n_low=nlo,n_high=nhi,
            density_jump=nhi-nlo,latent_energy=latent,
            fortran_maxwell_identity=identity,
            from_cubic_quintic_local_root_P=float(root[0]["p_t"]),
            local_root_pressure_difference=p-float(root[0]["p_t"]),
            raw_source_Crust_sha256=sha(ROOT/CAT[pair["crust_source_id"]]["compose_path"]/"eos.thermo"),
            raw_source_Core_sha256=sha(ROOT/CAT[pair["core_source_id"]]["compose_path"]/"eos.thermo"),
            interpolated_source_component_rows=[lc,hc],
            source_interpolation_diagnostics=[rl,rh],matched_rows=matched.n,
            serialized_readback_abs=maxdiff,
            output_hashes={z.name:sha(z) for z in target.glob("eos.*")},
            native_report=report)
    except Exception as e:
        result.update(status="BLOCKED_NATIVE_MAXWELL",
                      error_type=type(e).__name__,reason=str(e))
    save(dest/"RESULT.json",result)
    print("NATIVE_MAXWELL_2X2",tier,crust,core,order,result["status"],
          result.get("transition_P"),flush=True)
    return dict(tier=tier,crust=crust,core=core,order=order,status=result["status"])

def main():
    assert MATCHING_METHODS==("maxwell","css") and INTERPOLATION_METHODS==("cubic","quintic")
    catalog=json.loads((OUT/"CAMPAIGN.json").read_text())
    cases=[x for x in catalog["cases"] if x["both_status"]=="FINITE_ONLY_NOT_WHOLE"]
    if len(cases)!=7:raise ValueError("frozen seven native finite pairs required")
    jobs=[(pair,order) for pair in cases for order in INTERPOLATION_METHODS]
    if DEST.exists():raise FileExistsError("existing native Maxwell reference evidence")
    DEST.mkdir(parents=True)
    started=time.monotonic()
    with ThreadPoolExecutor(max_workers=4) as pool:
        res=[f.result() for f in as_completed([pool.submit(one,*args) for args in jobs])]
    rec=dict(status="COMPLETE_WITH_PRESERVED_REJECTIONS",matching="maxwell",
             interpolation_orders=list(INTERPOLATION_METHODS),jobs=len(jobs),
             parallel_workers=4,
             native_maxwell_accepted=sum(z["status"]=="PASS_NATIVE_MAXWELL_MATCHED" for z in res),
             duration_seconds=time.monotonic()-started,
             results=sorted(res,key=lambda z:(z["tier"],z["crust"],z["core"],z["order"])))
    save(DEST/"CAMPAIGN.json",rec)
    print("NATIVE_MAXWELL_COMPLETED",rec["native_maxwell_accepted"],"/",rec["jobs"])
if __name__=="__main__":main()
