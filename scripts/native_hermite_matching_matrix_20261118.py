"""Actual native cubic_F/quintic_F/cubic_P/quintic_P/CSS matching matrix.

Shared two crust and two core source inputs for each method. Maxwell root is
a *separate* construction and is never mislabeled CSS or smooth Hermite.
"""
from __future__ import annotations
import sys,json,hashlib,time,traceback
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos import read_compose,write_compose,match_crust_core
from compose_eos.native import executable
METHODS=("cubic_F","quintic_F","cubic_P","quintic_P","css")
PAIRS=(("SPG_M1","RG_SLY4"),("SPG_M3","RG_SLY4"),("SPG_M1","PK1_CANONICAL"),("SPG_M3","PK1_CANONICAL"))
OUT=ROOT/"data/cubic_quintic_matching_20261118/native_matching"
def source(role,name):
    if role=="crust":
        return ROOT/f"data/research_matching_20261117/spg/native_129_repaired/{name}/generated"
    if name=="RG_SLY4":return ROOT/"data/native_dense_interpolation_20261118/core/RG_SLY4/RG_SLY4"
    if name=="PK1_CANONICAL":return ROOT/"data/native_dense_interpolation_20261118/rmf/PK1_CANONICAL/table/compose"
    raise ValueError(name)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):
    p.parent.mkdir(parents=True,exist_ok=True)
    def conv(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        raise TypeError(str(type(v)))
    p.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False,default=conv)+"\n")
def run_pair(pair):
    crust,core=pair
    dcr=source("crust",crust);dco=source("core",core)
    folder=OUT/f"{crust}__{core}"
    if folder.exists():raise FileExistsError(folder)
    folder.mkdir(parents=True)
    l=read_compose(dcr);h=read_compose(dco)
    source_files={tag:{q.name:sha(q) for q in sorted(p.glob("eos.*"))} for tag,p in (("crust",dcr),("core",dco))}
    result=dict(pair=list(pair),input_hashes=source_files,raw_source_density_intervals={"crust":[float(l.nb[0]),float(l.nb[-1])],"core":[float(h.nb[0]),float(h.nb[-1])]},
      native_executable_sha256=sha(executable()),
      native_methods=list(METHODS),model_scope="matched finite-domain barotropic EOS from supplied native sources, not unique unified ground state",
      source_eos_file_present={"crust":l.composition_file_present,"core":h.composition_file_present},
      physical_energy_shift="NONE",status="RUNNING",methods={})
    save(folder/"RESULT.json",result)
    for m in METHODS:
        target=folder/m
        target.mkdir()
        t=time.monotonic()
        try:
            matched,info=match_crust_core(l,h,crust,core,m,40,drop_composition=True)
            write_compose(target/"compose",matched)
            back=read_compose(target/"compose")
            diff={k:float(np.max(np.abs(getattr(back,k)-getattr(matched,k)))) for k in ("nb","Q1","Q3","Q6","Q7")}
            if max(diff.values())>1e-11:raise ValueError("serialized matched CompOSE readback differs")
            if np.any(np.diff(back.nb)<=0):raise ValueError("matching output density order invalid")
            if not info["consistency"]["overall"]["pass_"]:raise ValueError("native match failed overall consistency")
            summary_keys=("i1","i2","n1","n2","P1","P2","mu1","mu2","eps1","eps2",
                "cs2_1","cs2_2","cs2_left_br","cs2_right_br","cs2_min_bridge","cs2_max_bridge",
                "first_law_relative_residual","is_exact_maxwell","css_single_sound_speed",
                "css_transition_nb_low","css_transition_nb_high","css_transition_delta_eps",
                "css_transition_pressure","css_transition_mu","css_transition_low_index",
                "css_transition_high_index","css_core_bridge_index",
                "phase_transition_class","construction_kind","regularity_class","seidov_ratio")
            snapshot={k:info.get(k) for k in summary_keys}
            rec=dict(status="CONSTRUCTED_VERIFIED_NATIVE",actual_native_method=m,
              n_input=[l.n,h.n],output_rows=back.n,
              transition_or_bridge=snapshot,serialized_readback_max_abs=diff,
              output_files_sha256={q.name:sha(q) for q in sorted((target/"compose").glob("eos.*"))},
              native_report=info,source_unchanged=True)
        except Exception as exc:
            rec=dict(status="REJECTED_OR_UNSUPPORTED",exception_type=type(exc).__name__,reason=str(exc),
                     traceback=traceback.format_exc()[-2500:])
        rec["elapsed_seconds"]=time.monotonic()-t
        result["methods"][m]=rec
        save(folder/"RESULT.json",result)
        print("MATCH",crust,core,m,rec["status"],str(rec.get("reason",""))[:180],flush=True)
    result["status"]="COMPLETE_WITH_PRESERVED_REJECTIONS"
    save(folder/"RESULT.json",result)
    return result
def main():
    from concurrent.futures import ThreadPoolExecutor,as_completed
    start=time.monotonic()
    outcomes=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for future in as_completed([pool.submit(run_pair,pair) for pair in PAIRS]):
            z=future.result();outcomes.append(dict(pair=z["pair"],
              generated=sum(x["status"]=="CONSTRUCTED_VERIFIED_NATIVE" for x in z["methods"].values()),
              rejections={k:v.get("reason") for k,v in z["methods"].items() if v["status"]!="CONSTRUCTED_VERIFIED_NATIVE"}))
    save(OUT/"CAMPAIGN.json",dict(status="EXECUTED",pairs=outcomes,elapsed_seconds=time.monotonic()-start,
        requested_combinations=len(PAIRS)*len(METHODS),completed=sum(q["generated"] for q in outcomes)))
    print("MATCHING_CAMPAIGN_DONE",len(outcomes),sum(q["generated"] for q in outcomes))
if __name__=="__main__":main()
