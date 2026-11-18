"""Compare EOS interpolation against *newly calculated* native dense states.

This is not a comparison against the true neutron-star EOS. It is a
source-fixed, phase-separated holdout benchmark against actual native
Fortran model evaluations. Physical sound speed uses finite-difference
reference with separate cross-derivative diagnostic.
"""
from __future__ import annotations
import argparse,json,hashlib,sys,traceback,time
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicHermiteSpline,PchipInterpolator,CubicSpline,Akima1DInterpolator
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos.api import read_compose
from compose_eos.whole_crust import CausalExponentialInterval
from matching_interpolation_shapes_20261118 import CausalShapeInterval,SlopeMatchedBernsteinInterval
OUT=ROOT/"data/native_dense_interpolation_20261118"
RMF_MODELS=("NL3_PROVIDER_PRECISION","PK1_CANONICAL","TM1_CANONICAL",
           "DDMEX_PROVIDER_ROUNDED_N0","DDMEX_ORIGINAL_TABLE_RHO0")
CORE_MODELS=("APR_APR","BBB_BHF2","GM_GM1","RG_KDE0v1","RG_SLY4")
GLOBAL_METHODS=("hermite","pchip","cubic","akima","linear")
SHAPE_METHODS=("exponential","rational","power","bernstein")
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    def conv(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        raise TypeError(str(type(v)))
    p.write_text(json.dumps(obj,indent=2,allow_nan=False,default=conv)+"\n")
def load_native(family,model):
    if family=="rmf":
        if model not in RMF_MODELS:raise ValueError("unknown RMF model")
        path=OUT/"rmf"/model/"table"/"compose"
        receipt=json.loads((OUT/"rmf"/model/"RESULT.json").read_text())
        assert receipt["status"]=="GENERATED_VERIFIED_DIRECT_NATIVE"
    elif family=="core":
        if model not in CORE_MODELS:raise ValueError("unknown core model")
        path=OUT/"core"/model/model
        receipt=json.loads((OUT/"core"/model/"SUMMARY.json").read_text())
        assert receipt["status"]=="COMPLETE"
    else:raise ValueError("unknown family")
    eos=read_compose(path)
    n=np.asarray(eos.nb,dtype=float);mu=np.asarray(eos.mub,dtype=float)
    p=np.asarray(eos.p,dtype=float);e=np.asarray(eos.eps,dtype=float)
    phase=np.asarray(eos.Iphase,dtype=int)
    if not(len(n)>3 and np.all(np.diff(n)>0) and np.all(np.diff(mu)>0)
            and np.all(np.diff(p)>0) and np.all(e>0) and np.all(mu>0)):
        raise ValueError("native source is not strictly monotone in n,mu,P")
    return dict(family=family,model=model,n=n,mu=mu,p=p,epsilon=e,
        phase=phase,sha256=sha(path/"eos.thermo"),
        composition_scope="beta_equilibrated native homogeneous/phase-specified output",
        source_table_dir=str(path.relative_to(ROOT)))
def phase_runs(phase):
    labels=np.asarray(phase)
    if labels.ndim!=1 or not len(labels):raise ValueError("one-dimensional nonempty phase labels required")
    cuts=np.r_[0,1+np.flatnonzero(np.diff(labels)!=0),len(labels)]
    return [(int(a),int(b)) for a,b in zip(cuts[:-1],cuts[1:])]
def train_holdout_indices(length,stride):
    if length<12 or stride<2:raise ValueError("need at least twelve source rows and stride>=2")
    tr=np.unique(np.r_[np.arange(0,length,stride),length-1])
    hold=np.setdiff1d(np.arange(length),tr,assume_unique=True)
    assert np.intersect1d(tr,hold).size==0
    return tr,hold
def predict_global(method,src,tr,hold):
    xx=src["mu"][tr];yy=src["p"][tr];xxh=src["mu"][hold]
    if method=="linear":
        slot=np.clip(np.searchsorted(xx,xxh,side="right")-1,0,len(xx)-2)
        slope=np.diff(yy)/np.diff(xx)
        pred_p=yy[slot]+slope[slot]*(xxh-xx[slot]);pred_n=slope[slot]
        dn=np.zeros_like(pred_n);pred_cs=np.full_like(pred_n,np.nan)
    else:
        if method=="hermite":spl=CubicHermiteSpline(xx,yy,src["n"][tr])
        elif method=="pchip":spl=PchipInterpolator(xx,yy)
        elif method=="cubic":spl=CubicSpline(xx,yy,bc_type="natural")
        elif method=="akima":spl=Akima1DInterpolator(xx,yy,method="akima")
        else:raise ValueError(method)
        pred_p=spl(xxh);pred_n=spl(xxh,nu=1);dn=spl(xxh,nu=2)
        with np.errstate(divide="ignore",invalid="ignore"):
            pred_cs=pred_n/(xxh*dn)
    return dict(mu=xxh,p=np.asarray(pred_p),n=np.asarray(pred_n),
        epsilon=xxh*np.asarray(pred_n)-np.asarray(pred_p),
        cs2=np.asarray(pred_cs),dn_dmu=np.asarray(dn))
def prediction_from_shape(shape,source,train,hold):
    fields=("p","n","epsilon","cs2","dn_dmu")
    pred={k:np.full(len(hold),np.nan) for k in fields}
    pred["mu"]=source["mu"][hold]
    original_index={int(h):i for i,h in enumerate(hold)}
    statistics=dict(segments_considered=0,segments_admitted=0,segments_rejected=0,
                    rejected_by_reason={},failure_examples=[])
    xx=source["mu"]
    if shape=="bernstein":
        with np.errstate(divide="ignore",invalid="ignore"):
            deriv=np.gradient(source["n"][train],xx[train],edge_order=2)
            cstr=source["n"][train]/(xx[train]*deriv)
    for i in range(len(train)-1):
        left,right=train[i:i+2]
        ids=[int(j) for j in hold if left<j<right]
        if not ids:continue
        statistics["segments_considered"]+=1
        values=[source[k][j] for j in (left,right) for k in ("n","p","epsilon")]
        try:
            if shape=="exponential":
                b=CausalExponentialInterval(*values)
            elif shape in ("power","rational"):
                b=CausalShapeInterval(*values,family=shape)
            elif shape=="bernstein":
                if (not np.isfinite(cstr[i]) or not np.isfinite(cstr[i+1])
                    or not 0<cstr[i]<=1 or not 0<cstr[i+1]<=1):
                    raise ValueError("noncausal_or_unresolved_training_slope")
                b=SlopeMatchedBernsteinInterval(*values,left_cs2=float(cstr[i]),
                  right_cs2=float(cstr[i+1]))
            else:raise ValueError("unknown shape")
            verdict=b.validate()
            if verdict["status"]!="PASS":
                raise ValueError("source_fixed_causal_interval_gate_failed")
            observed=b.evaluate_mu(xx[ids])
            if not np.all(np.isfinite(np.asarray(observed["epsilon"]))):
                raise ValueError("nonfinite interpolation")
            for name in ("p","n","epsilon"):
                for j,value in zip(ids,np.atleast_1d(observed[name])):
                    pred[name][original_index[j]]=value
            for j,value in zip(ids,np.atleast_1d(observed["cs2_over_c2"])):
                pred["cs2"][original_index[j]]=value
                pred["dn_dmu"][original_index[j]]=observed["n"][list(ids).index(j)]/(xx[j]*value)
            statistics["segments_admitted"]+=1
        except (ValueError,OverflowError,ZeroDivisionError,FloatingPointError) as exc:
            reason=str(exc)[:160]
            statistics["segments_rejected"]+=1
            statistics["rejected_by_reason"][reason]=statistics["rejected_by_reason"].get(reason,0)+1
            if len(statistics["failure_examples"])<6:
                statistics["failure_examples"].append(dict(lower_n=float(source["n"][left]),upper_n=float(source["n"][right]),reason=reason))
    return pred,statistics
def native_truth_cs2(s):
    n,mu,p,e=[s[k] for k in ("n","mu","p","epsilon")]
    dp=np.gradient(p,n,edge_order=2)
    de=np.gradient(e,n,edge_order=2)
    dn=np.gradient(n,mu,edge_order=2)
    with np.errstate(divide="ignore",invalid="ignore"):
        a=dp/de;b=n/(mu*dn)
    return dict(cs2_from_P_over_e=a,cs2_from_n_over_mu=b,
                max_abs_cross_form=float(np.nanmax(np.abs(a-b))),
                frac_causal=float(np.mean((a>0)&(a<=1))),
                count_acausal=int(np.sum(a>1)),
                acausal_density_fm3=[float(n[i]) for i in np.flatnonzero(a>1)[[0,-1]]] if np.any(a>1) else [])
def summarize(method,source,hold,pred,truth,meta=None):
    p=np.asarray(pred["p"]);n=np.asarray(pred["n"]);e=np.asarray(pred["epsilon"]);mu=source["mu"][hold]
    c=np.asarray(pred["cs2"]);nc=truth["cs2_from_P_over_e"][hold]
    native_valid=np.isfinite(nc)&(nc>0)&(nc<=1)
    valid=np.isfinite(p)&np.isfinite(n)&np.isfinite(e)&(n>0)
    area=valid & native_valid
    cs_good=area & np.isfinite(c)&(c>0)&(c<=1+5e-12)
    def stats(values,mask):
        v=np.abs(np.asarray(values)[mask])
        if not len(v):return dict(count=0,rms=None,p95=None,max=None)
        return dict(count=int(len(v)),rms=float(np.sqrt(np.mean(v*v))),
                    p95=float(np.quantile(v,.95)),max=float(np.max(v)))
    prec=(p-source["p"][hold])/np.maximum(np.abs(source["p"][hold]),1e-3)
    nrec=(n-source["n"][hold])/source["n"][hold]
    edev=(e-source["epsilon"][hold])/source["n"][hold]
    csdev=c-nc
    res=dict(holdout_count=int(len(hold)),native_causal_holdout=int(np.sum(native_valid)),
         finite_prediction_fraction=float(np.mean(valid)),physical_among_native_causal=float(np.sum(cs_good)/max(np.sum(native_valid),1)),
         invalid_derived_cs2_count=int(np.sum(area & ~(np.isfinite(c)&(c>0)&(c<=1+5e-12)))),
         max_first_law_Euler_abs=float(np.max(abs(mu[valid]*n[valid]-p[valid]-e[valid]))) if np.any(valid) else None,
         pressure_relative=stats(prec,valid),number_density_relative=stats(nrec,valid),
         energy_per_baryon_MeV=stats(edev,valid),sound_speed_squared_absolute=stats(csdev,cs_good))
    if meta is not None:res["interpolation_diagnostics"]=meta
    return res
def generate_one(family,model):
    start=time.monotonic();s=load_native(family,model)
    resultdir=OUT/"cv"/family/model
    if resultdir.exists():raise FileExistsError("CV output would overwrite existing evidence")
    resultdir.mkdir(parents=True)
    result=dict(schema="ns_cowling.native_dense_interpolation_holdout.v1",
       family=family,model=model,status="RUNNING",
       native_source_sha256=s["sha256"],native_source_table=s["source_table_dir"],
       source_n_nodes=int(len(s["n"])),n_range_fm3=[float(s["n"][0]),float(s["n"][-1])],
       source_phase_ids=sorted(set(map(int,s["phase"]))),
       protocol="phase-separated withheld native evaluation points; no random splitting; no reinterpretation of F/epsilon",
       comparison_coordinate="P(mu), n(mu), epsilon(mu) at native test chemical potentials; logP(n) controls not pooled",
       numerical_causality="native direct finite differences distinguish original source acausality from interpolation violations",
       test_density_strides=[16,8] if family=="rmf" else [8,4],runs=[])
    save(resultdir/"RESULT.json",result)
    phase=phase_runs(s["phase"])
    for a,b in phase:
        branch={k:v[a:b] if isinstance(v,np.ndarray) else v for k,v in s.items()}
        if b-a<16:
            result["runs"].append(dict(status="SKIPPED_SHORT_NATIVE_PHASE",phase=int(s["phase"][a]),rows=b-a))
            continue
        direct=native_truth_cs2(branch)
        for stride in result["test_density_strides"]:
            tr,hold=train_holdout_indices(len(branch["n"]),stride)
            raw=dict(phase=int(branch["phase"][0]),source_row_range=[a,b],train_stride=stride,
                training_size=len(tr),heldout_size=len(hold),native_physics=direct,
                metrics={},status="RUNNING")
            result["runs"].append(raw)
            save(resultdir/"RESULT.json",result)
            for method in GLOBAL_METHODS+SHAPE_METHODS:
                try:
                    if method in GLOBAL_METHODS:
                        predicted=predict_global(method,branch,tr,hold);extra=None
                    else:
                        predicted,extra=prediction_from_shape(method,branch,tr,hold)
                    metrics=summarize(method,branch,hold,predicted,direct,extra)
                    metrics["status"]="COMPUTED_WITH_NATIVE_HOLDOUT"
                    raw["metrics"][method]=metrics
                    print("BENCH",family,model,"phase",raw["phase"],"stride",stride,
                          method,"valid",f"{metrics['finite_prediction_fraction']:.3f}",
                          "physical",f"{metrics['physical_among_native_causal']:.3f}",
                          "n_rms",metrics["number_density_relative"]["rms"],flush=True)
                except Exception as exc:
                    raw["metrics"][method]=dict(status="BLOCKED_METHOD",type=type(exc).__name__,
                                                message=str(exc),traceback=traceback.format_exc()[-1800:])
                    print("BENCH_FAILED",family,model,method,type(exc).__name__,str(exc)[:110],flush=True)
                save(resultdir/"RESULT.json",result)
            raw["status"]="COMPLETE_WITH_REJECTED_METHODS"
            save(resultdir/"RESULT.json",result)
    result["duration_seconds"]=time.monotonic()-start
    result["status"]="COMPLETED_SCOPED_WITH_NATIVE_HOLDOUT"
    save(resultdir/"RESULT.json",result)
    print("BENCH_COMPLETED",family,model,len(result["runs"]),f"{result['duration_seconds']:.1f}",flush=True)
    return result
if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--family",choices=("rmf","core"),required=True)
    parser.add_argument("--model",required=True)
    args=parser.parse_args()
    generate_one(args.family,args.model)
