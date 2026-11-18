"""Independently re-search Maxwell pressure crossings from untouched native source pairs.

Cubic/quintic data fits within each original phase, never across an existing
Maxwell jump. Coarse fits are compared to the native full-source Maxwell root.
CSS is a distinct construction; comparison to Maxwell is a model-choice
difference, NOT a Hermite interpolation error.
"""
from pathlib import Path
import sys,hashlib,json,time,traceback
import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import CubicHermiteSpline,PchipInterpolator
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from native_hermite_matching_matrix_20261118 import PAIRS,source
from compare_hermite_orders_20261118 import quintic_hermite
from compose_eos import read_compose
from compose_eos.api import construct_maxwell_transition,_thermo_only_view
OUT=ROOT/"data/cubic_quintic_matching_20261118/maxwell"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    def conv(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        raise TypeError(str(type(v)))
    p.write_text(json.dumps(x,indent=2,allow_nan=False,default=conv)+"\n")
def family_fit(d,method,stride):
    mu=np.asarray(d.mub,float);p=np.asarray(d.p,float);n=np.asarray(d.nb,float)
    tr=np.unique(np.r_[np.arange(0,len(mu),stride),len(mu)-1])
    x,y,z=mu[tr],p[tr],n[tr]
    if np.any(np.diff(x)<=0) or np.any(np.diff(z)<=0):raise ValueError("source component monotonicity fails")
    if method=="cubic":
        spline=CubicHermiteSpline(x,y,z)
    elif method=="quintic":
        susceptibility=PchipInterpolator(x,z).derivative()(x)
        spline=quintic_hermite(x,y,z,susceptibility)
    else:raise ValueError(method)
    return dict(spline=spline,mu=[float(x[0]),float(x[-1])],retained=int(len(x)),stride=stride)
def scan_one(crust,core,method,stride_low,stride_high):
    l=family_fit(crust,method,stride_low)
    h=family_fit(core,method,stride_high)
    lower=max(l["mu"][0],h["mu"][0]);upper=min(l["mu"][1],h["mu"][1])
    if not upper>lower:raise ValueError("no source overlap in Euler chemical potential")
    # independent fixed scan. Piecewise polynomials can produce multiple crossings.
    mu=np.linspace(lower,upper,3001)
    residual=h["spline"](mu)-l["spline"](mu)
    brackets=[(float(mu[i]),float(mu[i+1])) for i in range(len(mu)-1) if residual[i]*residual[i+1]<0]
    roots=[]
    for a,b in brackets:
        root=float(brentq(lambda m:float(h["spline"](m)-l["spline"](m)),a,b,xtol=1e-11))
        nlo=float(l["spline"](root,nu=1));nhi=float(h["spline"](root,nu=1))
        pressure=float(.5*(l["spline"](root)+h["spline"](root)))
        delta_eps=root*(nhi-nlo)
        speed=[]
        for s,n in ((l["spline"],nlo),(h["spline"],nhi)):
            chi=float(s(root,nu=2))
            speed.append(float(n/(root*chi)) if chi!=0 else None)
        roots.append(dict(mu_MeV=root,P_t_MeV_fm3=pressure,
                          n_low=nlo,n_high=nhi,density_jump=nhi-nlo,
                          delta_epsilon_MeV_fm3=delta_eps,
                          c_s2_one_sided=speed,
                          positive_density_jump=nhi>nlo,
                          lower_pressure_dominant_before=bool(float(l["spline"](a))>float(h["spline"](a))),
                          higher_pressure_dominant_after=bool(float(h["spline"](b))>float(l["spline"](b))),
                          local_pressure_residual=float(h["spline"](root)-l["spline"](root))))
    return dict(method=method,crust_training_stride=stride_low,core_training_stride=stride_high,
                crust_training_rows=l["retained"],core_training_rows=h["retained"],
                common_mu_MeV=[lower,upper],roots=roots,crossings=len(roots),
                globally_accepted_single_root=len(roots)==1 and roots[0]["positive_density_jump"]
                   and roots[0]["lower_pressure_dominant_before"] and roots[0]["higher_pressure_dominant_after"])
def run_pair(pair):
    c,k=pair
    cr=read_compose(source("crust",c));co=read_compose(source("core",k))
    folder=OUT/f"{c}__{k}"
    if folder.exists():raise FileExistsError(folder)
    folder.mkdir(parents=True)
    result=dict(pair=[c,k],status="RUNNING",physical_source="unmodified two native CompOSE models; shared physical reference",
      source_hashes={side:{p.name:sha(p) for p in src.glob("eos.*")} for side,src in (("crust",source("crust",c)),("core",source("core",k)))},fits={})
    save(folder/"RESULT.json",result)
    try:
        native,rep=construct_maxwell_transition(_thermo_only_view(cr),_thermo_only_view(co),
           assume_local_charge_neutrality=True,return_report=True)
        sel={k:float(rep[k]) for k in ("mu_transition","pressure_transition","n_low","n_high","density_jump","latent_energy_density","maxwell_identity_residual")}
        result["native_full_grid_maxwell"]=dict(status="PASS_SCOPED",**sel,
                 output_rows=native.n,
                 native_success_not_global_geometry_proof=True)
    except Exception as e:
        result["native_full_grid_maxwell"]=dict(status="REJECTED",error=str(e),
                                                category=getattr(e,"category",None))
    for method in ("cubic","quintic"):
        for level,cs,ks in (("native_full",(1),(1)),("sparse_4x",(4),(8)),("sparse_8x",(8),(16))):
            key=f"{method}_{level}"
            try:
                fit=scan_one(cr,co,method,cs,ks)
                native=result["native_full_grid_maxwell"]
                if native["status"]=="PASS_SCOPED" and fit["globally_accepted_single_root"]:
                    a=fit["roots"][0];fit["comparison_native_maxwell"]=dict(
                      mu_relative=(a["mu_MeV"]-native["mu_transition"])/native["mu_transition"],
                      pressure_relative=(a["P_t_MeV_fm3"]-native["pressure_transition"])/native["pressure_transition"],
                      density_jump_relative=(a["density_jump"]-native["density_jump"])/native["density_jump"],
                      latent_heat_relative=(a["delta_epsilon_MeV_fm3"]-native["latent_energy_density"])/native["latent_energy_density"])
                fit["status"]="ROOT_SCAN_COMPUTED"
            except Exception as e:
                fit=dict(status="REJECTED",error=str(e),category=type(e).__name__)
            result["fits"][key]=fit
            save(folder/"RESULT.json",result)
            print("MAXWELL",c,k,key,fit["status"],fit.get("crossings","?"),
                  fit.get("globally_accepted_single_root","?"),flush=True)
    result["status"]="COMPLETED_WITH_PRESERVED_REJECTIONS"
    save(folder/"RESULT.json",result)
    return result
def main():
    from concurrent.futures import ThreadPoolExecutor,as_completed
    OUT.mkdir(parents=True,exist_ok=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        results=[]
        for f in as_completed([pool.submit(run_pair,p) for p in PAIRS]):
            x=f.result()
            results.append(dict(pair=x["pair"],native_status=x["native_full_grid_maxwell"]["status"],
                  fits={k:v["status"] for k,v in x["fits"].items()}))
    save(OUT/"CAMPAIGN.json",dict(status="EXECUTED",results=results))
    print("ALL_MAXWELL_DONE",len(results))
if __name__=="__main__":main()
