"""Exactly two matching constructions with two interpolation orders.

Matching is ONLY:
    (i) true Maxwell coexistence P_-(mu)=P_+(mu), or
    (ii) the native CSS phase-transition matching contract.

Cubic and quintic are *input EOS representations*, never matching methods.
Native source rows are retained exactly and only strictly causal, same-phase
midpoints are generated. Source/provenance tier is unchanged; derived
representations are clearly labeled, not native-original thermodynamic data.
"""
from __future__ import annotations
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import json,hashlib,time,sys,traceback
import numpy as np
from scipy.interpolate import CubicHermiteSpline,PchipInterpolator,BPoly

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from compose_eos.eos_data import EOSData
from compose_eos import read_compose,match_crust_core,write_compose
from all_origin_maxwell_css_20261118 import OUT,sha
from check_hermite_roots_for_dual_pairs_20261118 import source,CAT

MATCHING_METHODS=("maxwell","css")
INTERPOLATION_METHODS=("cubic","quintic")
EXPECTED_VARIANTS=tuple(f"{method}+{matching}" for matching in MATCHING_METHODS for method in INTERPOLATION_METHODS)
OUTPUT=OUT/"HERMITED_INPUT_2X2"

def save(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    def conv(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    p.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False,default=conv)+"\n")


def _thermo_from_eos(template,n,p,e,mu,phase):
    n,p,e,mu=(np.asarray(z,float) for z in (n,p,e,mu))
    if not len(n)==len(p)==len(e)==len(mu) or len(n)<3:
        raise ValueError("no valid source points")
    if any(np.any(~np.isfinite(x)) or np.any(x<=0) for x in (n,p,e,mu)):
        raise ValueError("positive finite thermodynamic states required")
    if np.any(np.diff(n)<=0):
        raise ValueError("new source not density-ordered")
    z=np.zeros(len(n),dtype=float)
    return EOSData(nb=n,mn=template.mn,mp=template.mp,
        Q1=p/n,Q2=z.copy(),Q3=mu/template.mn-1.,Q4=z.copy(),Q5=z.copy(),
        Q6=e/(n*template.mn)-1.,Q7=e/(n*template.mn)-1.,
        q_extra=z.copy(),Iphase=np.asarray(phase,dtype=np.int64),
        Yn=z.copy(),Yp=z.copy(),Ye=z.copy(),Ymu=z.copy(),Aav=z.copy(),Zav=z.copy(),Yh=z.copy(),
        Il=template.Il,iT=template.iT,iYq=template.iYq,
        temperature_mev=template.temperature_mev,
        charge_fraction_grid=template.charge_fraction_grid,
        composition_file_present=False,composition_available=False,beta_equilibrated=template.beta_equilibrated,
        source_dir=None)


def resample_eos(source_data:EOSData,order:str):
    if order not in INTERPOLATION_METHODS:
        raise ValueError("only cubic/quintic interpolation orders")
    n=np.asarray(source_data.nb,float)
    p=np.asarray(source_data.p,float)
    e=np.asarray(source_data.eps,float)
    mu=np.asarray(source_data.mub,float)
    phase=np.asarray(source_data.Iphase,dtype=int)
    points=[];inserted=0;rejected=0;skipped_phase=0
    length=len(n)
    if length<3:raise ValueError("insufficient source knots")
    # Build bounded smooth source subsegments only. No cross-iphase blending and
    # no accepted-state stitching across a large unsampled density gap.
    cuts=np.r_[0,np.flatnonzero(np.diff(phase)!=0)+1,length]
    for a,b in zip(cuts[:-1],cuts[1:]):
        valid=np.array([(mu[i+1]>mu[i]>0 and n[i+1]>n[i]>0 and
                         p[i+1]>p[i]>0 and e[i+1]>e[i]>0 and
                         n[i+1]/n[i]<1.4)
                        for i in range(a,b-1)],bool)
        # Enumerate only maximal, monotone segments of >=4 source knots.
        segments=[];start=None
        for i,okay in enumerate(valid):
            if okay and start is None:start=a+i
            if not okay and start is not None:
                if a+i+1-start>=4:segments.append((start,a+i+1))
                start=None
        if start is not None and b-start>=4:segments.append((start,b))
        cache={}
        for left,right in segments:
            x=mu[left:right];y=p[left:right];density=n[left:right]
            if order=="cubic":
                spl=CubicHermiteSpline(x,y,density)
            else:
                chi=PchipInterpolator(x,density).derivative()(x)
                if not np.all(np.isfinite(chi)) or np.any(chi<=0):
                    continue
                spl=BPoly.from_derivatives(x,[[y[i],density[i],chi[i]] for i in range(len(x))])
            for i in range(left,right-1):
                umid=.5*(mu[i]+mu[i+1])
                pmid=float(spl(umid));nmid=float(spl(umid,nu=1));chi=float(spl(umid,nu=2))
                emid=umid*nmid-pmid
                cs2=nmid/(umid*chi) if chi!=0 else np.inf
                if (p[i]<pmid<p[i+1] and n[i]<nmid<n[i+1]
                    and e[i]<emid<e[i+1] and 0<cs2<=1 and np.isfinite(cs2)):
                    cache[i]=(nmid,pmid,emid,umid,int(phase[i]))
                else:rejected+=1
        for i in range(a,b):
            points.append((n[i],p[i],e[i],mu[i],int(phase[i])))
            if i in cache:
                points.append(cache[i]);inserted+=1
        skipped_phase+=max(0,b-a-1-len(cache))
    arr=np.asarray(points,float)
    if not np.all(np.diff(arr[:,0])>0):
        raise ValueError("phase-preserving inserted nodes broke density ordering")
    result=_thermo_from_eos(source_data,arr[:,0],arr[:,1],arr[:,2],arr[:,3],arr[:,4])
    receipt=dict(original_rows=len(n),interpolated_rows=len(points),inserted_rows=inserted,
                 rejected_intermediate_rows=rejected,unfilled_source_phase_intervals=skipped_phase,
                 method=order,matches_only_original_source_rows=True,
                 derived_input_representation=True,not_provider_original_rows=True,
                 physical_law="P(mu), n=dP/dmu, epsilon=mu*n-P at inserted rows",
                 physical_scope="strict same-phase causal interpolation at new midpoints only")
    return result,receipt


def run_one(entry):
    tier,crust,core=entry["tier"],entry["crust"],entry["core"]
    case=OUTPUT/tier/(crust+"__"+core)
    if case.exists():raise FileExistsError(case)
    case.mkdir(parents=True)
    pair=(tier,crust,core)
    known=next(x for x in json.loads((OUT/"HERMITE_ROOT_CROSSCHECK.json").read_text())["cases"]
               if (x["tier"],x["crust"],x["core"])==pair)
    native_css=entry["css"]
    original_l=source(entry["crust_source_id"])
    original_h=source(entry["core_source_id"])
    record=dict(tier=tier,crust=crust,core=core,
                original_source_ids=[entry["crust_source_id"],entry["core_source_id"]],
                source_nuclear_functionals_different=not entry["same_functional"],
                source_origin_never_mixed=True,methods=list(MATCHING_METHODS),
                interpolations=list(INTERPOLATION_METHODS),
                source_matched_native_CSS=entry["css"],
                original_legacy_maxwell=entry["maxwell"],
                interpolation_pair_results={},status="RUNNING")
    save(case/"RESULT.json",record)
    for order in INTERPOLATION_METHODS:
        # The two matching algorithms consume the very same input-order
        # representation and the output is named ONLY by matching type.
        maxwell=known["methods"][order]
        roots=maxwell["roots"]
        record["interpolation_pair_results"][order+"+maxwell"]=dict(
            status="PASS_SCOPED_LOCAL_MAXWELL_ROOT" if len(roots)==1 else "UNRESOLVED",
            input_interpolation_order=order,matching="maxwell",
            roots=roots,limitation="local source root; independent native Fortran Maxwell remains admission authority")
        try:
            lo,ra=resample_eos(original_l,order)
            hi,rb=resample_eos(original_h,order)
            if min(ra["inserted_rows"],rb["inserted_rows"])<=0:
                raise ValueError("not enough admitted intermediate source points for input interpolation test")
            matched,report=match_crust_core(lo,hi,crust,core,"css",40,drop_composition=True)
            if not report["consistency"]["overall"]["pass_"]:
                raise ValueError("derived representation fails native CSS consistency")
            target=case/order/"css"
            target.mkdir(parents=True)
            write_compose(target/"compose",matched)
            back=read_compose(target/"compose")
            assert matched.n==back.n
            delta=max(float(np.max(abs(getattr(back,k)-getattr(matched,k)))) for k in ("nb","Q1","Q3","Q6","Q7"))
            if delta>1e-9:
                raise ValueError(f"derived CSS serialisation changed: {delta}")
            css=dict(status="PASS_DERIVED_INTERPOLATED_SOURCE",
                matching="css",input_interpolation_order=order,source_crust_interpolant=ra,
                source_core_interpolant=rb,
                P_t=float(report["css_transition_pressure"]),
                density_jump=float(report["css_transition_nb_high"]-report["css_transition_nb_low"]),
                delta_epsilon=float(report["css_transition_delta_eps"]),
                cs2_css=float(report["css_single_sound_speed"]),
                output_n=int(back.n),source_remained_immutable=True,
                source_provenance_tier=tier,
                source_representation="thermodynamically derived interpolation from same-tier original knots; NOT raw original provider numerical data",
                actual_native_matching="css",
                serialized_readback_abs=delta)
        except Exception as exc:
            css=dict(status="BLOCKED_OR_REJECTED",matching="css",input_interpolation_order=order,
                     error_type=type(exc).__name__,error=str(exc),
                     source_provenance_tier=tier,
                     source_representation="in-memory derived interpolated EOS; not original")
        record["interpolation_pair_results"][order+"+css"]=css
        save(case/"RESULT.json",record)
        print("INTERPOLATION_2X2",tier,crust,core,order,
              "maxwell",record["interpolation_pair_results"][order+"+maxwell"]["status"],
              "css",css["status"],flush=True)
    record["status"]="COMPLETED_WITH_PRESERVED_REJECTIONS"
    save(case/"RESULT.json",record)
    return record


def main():
    archive=json.loads((OUT/"CAMPAIGN.json").read_text())
    inputs=[x for x in archive["cases"] if x["both_status"]=="FINITE_ONLY_NOT_WHOLE"]
    if len(inputs)!=7:raise ValueError("selected seven native finite pairs required")
    OUTPUT.mkdir(parents=True,exist_ok=True)
    with ThreadPoolExecutor(max_workers=3) as pool:
        result=[t.result() for t in as_completed([pool.submit(run_one,x) for x in inputs])]
    summary=dict(status="FINISHED_SCOPED",matchings=list(MATCHING_METHODS),
                interpolation_orders=list(INTERPOLATION_METHODS),
                inputs_count=len(inputs),variants_expected=7*4,
                outcomes=[dict(tier=x["tier"],crust=x["crust"],core=x["core"],
                    variants={k:v["status"] for k,v in x["interpolation_pair_results"].items()}) for x in result])
    save(OUTPUT/"CAMPAIGN.json",summary)
    print("INTERPOLATION_2X2_DONE",len(result),
          "VARIANTS",sum(len(x["interpolation_pair_results"]) for x in result))
if __name__=="__main__":main()
