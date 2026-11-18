"""Independent cubic/quintic P(mu) root consistency on the seven dual finite pairs.

Source-retained (P,mu,n) are never modified; the quintic chi is evaluated only
from the SAME source component's n(mu) knots. No interpolation across phase
component gaps. This does not turn the CSS construction into a Hermite method.
"""
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import CubicHermiteSpline,PchipInterpolator,BPoly
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from compose_eos import read_compose
from compose_eos.maxwell_v2 import _segment_eosdata
from all_origin_maxwell_css_20261118 import OUT,catalog,bracket_roots
CAT={z["source_id"]:z for z in catalog()["sources"]}
DATA=json.loads((OUT/"CAMPAIGN.json").read_text())
SUCCESS=[z for z in DATA["cases"] if z["both_status"]=="FINITE_ONLY_NOT_WHOLE"]
LEGACY=json.loads((ROOT/"data/maxwell_three_families_20261117/campaign_v2/CAMPAIGN.json").read_text())
def source(srcid):
    z=CAT[srcid]
    d=read_compose(ROOT/z["compose_path"])
    if z["provenance_class"]=="new_public":
        a=z["selection"];d=_segment_eosdata(d,(0,int(a["source_start_zero_based"]),int(a["source_stop_exclusive"])))
    return d
def comps(d,entry,crust):
    rows=[]
    if entry.get("components"):
        for c in entry["components"]:
            a=int(c["start_zero_based"]);b=int(c["stop_exclusive"])
            if b-a>=4:rows.append((a,b,c["id"]))
    else:
        a=0
        # Source-preliminary phase distinctions never cross a phase label.
        phase=np.asarray(d.Iphase)
        end=np.r_[np.flatnonzero(np.diff(phase)!=0)+1,len(phase)]
        for b in end:
            if b-a>=4:rows.append((int(a),int(b),f"full_phase_{a}_{b}"))
            a=b
    return rows
def fit(d,sl,method):
    a,b,label=sl
    x=np.asarray(d.mub[a:b],float);p=np.asarray(d.p[a:b],float);n=np.asarray(d.nb[a:b],float)
    if np.any(np.diff(x)<=0) or np.any(np.diff(n)<=0):
        raise ValueError("source monotonicity")
    if method=="cubic":
        return CubicHermiteSpline(x,p,n)
    chi=PchipInterpolator(x,n).derivative()(x)
    if np.any(chi<=0):
        raise ValueError("nonpositive quintic source second derivatives")
    return BPoly.from_derivatives(x,[[p[i],n[i],chi[i]] for i in range(len(x))])
def root_scan(l,h,a,b,method,stride=1):
    pl=fit(l,a,method);ph=fit(h,b,method)
    lo=max(float(l.mub[a[0]]),float(h.mub[b[0]]))
    hi=min(float(l.mub[a[1]-1]),float(h.mub[b[1]-1]))
    if hi<=lo:return []
    roots=bracket_roots(pl,ph,lo,hi,1001)
    result=[]
    for mu in roots:
        p1,p2=float(pl(mu)),float(ph(mu))
        n1,n2=float(pl(mu,nu=1)),float(ph(mu,nu=1))
        if min(n1,n2)<=0 or n2<=n1:continue
        chi1,chi2=float(pl(mu,nu=2)),float(ph(mu,nu=2))
        if min(chi1,chi2)<=0:continue
        c1=n1/(mu*chi1);c2=n2/(mu*chi2)
        if not (0<c1<=1 and 0<c2<=1):continue
        dt=min(1e-4,1e-3*(hi-lo))
        if mu-dt<=lo or mu+dt>=hi:continue
        if not(float(pl(mu-dt))>float(ph(mu-dt)) and float(ph(mu+dt))>float(pl(mu+dt))):continue
        result.append(dict(mu_t=mu,p_t=(p1+p2)/2,n_low=n1,n_high=n2,
                           density_jump=n2-n1,latent_heat=mu*(n2-n1),
                           c2_sides=[c1,c2],
                           source_components=[a[2],b[2]],
                           P_abs_residual=abs(p1-p2)))
    return result
def main():
    rows=[]
    for z in SUCCESS:
        l=source(z["crust_source_id"]);h=source(z["core_source_id"])
        al=comps(l,CAT[z["crust_source_id"]],True)
        ah=comps(h,CAT[z["core_source_id"]],False)
        res=dict(tier=z["tier"],crust=z["crust"],core=z["core"],status="RUNNING",
                 same_source_exact_original=True,methods={})
        for method in ("cubic","quintic"):
            successes=[];rejections=[]
            for a in al:
                for b in ah:
                    if l.mub[a[1]-1]<h.mub[b[0]] or h.mub[b[1]-1]<l.mub[a[0]]:
                        continue
                    try:successes.extend(root_scan(l,h,a,b,method))
                    except Exception as e:
                        rejections.append(dict(components=[a[2],b[2]],why=str(e)[:120]))
            res["methods"][method]=dict(status="POSITIVE_ROOTS_FOUND" if successes else "NO_APPROVED_SPARSE_ROOT",
                                        roots=successes,components_rejected=rejections)
        res["status"]="COMPLETE"
        rows.append(res)
        print("HERMITE_7",res["tier"],res["crust"],res["core"],
              [(key,len(val["roots"])) for key,val in res["methods"].items()],flush=True)
    path=OUT/"HERMITE_ROOT_CROSSCHECK.json"
    path.write_text(json.dumps(dict(status="COMPLETE",eligible_dual_pair_count=len(SUCCESS),cases=rows),
                               indent=2,allow_nan=False)+"\n")
if __name__=="__main__":main()
