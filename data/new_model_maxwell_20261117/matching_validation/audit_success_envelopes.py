"""Independent polynomial check of every successful finite-component envelope."""
from pathlib import Path
import json
import sys
import numpy as np
from numpy.polynomial import Polynomial
from scipy.interpolate import CubicHermiteSpline

ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from maxwell_three_family_campaign import load, save, read_prepared, sha, rel
from compose_eos.maxwell_v2 import _segment_eosdata
BASE=ROOT/"data/new_model_maxwell_20261117"


def interval_poly(spline,a,b):
    i=min(max(np.searchsorted(spline.x,a,side="right")-1,0),len(spline.x)-2)
    c=spline.c[:,i]
    x=Polynomial([a-spline.x[i],b-a])
    return ((c[0]*x+c[1])*x+c[2])*x+c[3]


def unit_roots(poly):
    return sorted(float(z.real) for z in poly.roots()
        if abs(z.imag)<=1e-9 and -1e-10<=z.real<=1+1e-10)


def audit(low,high,target):
    ml=(low.eps+low.p)/low.nb
    mh=(high.eps+high.p)/high.nb
    spl=CubicHermiteSpline(ml,low.p,low.nb)
    sph=CubicHermiteSpline(mh,high.p,high.nb)
    lo,hi=max(ml[0],mh[0]),min(ml[-1],mh[-1])
    union=np.unique(np.r_[lo,ml[(ml>lo)&(ml<hi)],mh[(mh>lo)&(mh<hi)],hi])
    roots=[]
    stationary=0
    for a,b in zip(union[:-1],union[1:]):
        p=interval_poly(sph,a,b)-interval_poly(spl,a,b)
        stationary+=len(unit_roots(p.deriv()))
        for r in unit_roots(p):
            root=float(a+(b-a)*r)
            if not roots or abs(root-roots[-1])>1e-8:roots.append(root)
    first=float(sph(lo)-spl(lo))
    last=float(sph(hi)-spl(hi))
    mu=target["mu_transition"]
    slope=float(sph(mu,1)-spl(mu,1))
    error=abs(roots[0]-mu) if len(roots)==1 else None
    pt=target["pressure_transition"]
    residual=float(abs(sph(mu)-spl(mu)))
    passed=(len(roots)==1 and first<0 and last>0 and slope>0 and
        error <= 1e-10*max(1,abs(mu)) and residual<=1e-10*max(abs(pt),1e-300))
    return {"status":"PASS" if passed else "FAIL","common_mu_MeV":[float(lo),float(hi)],
        "union_intervals":len(union)-1,"stationary_points":stationary,
        "independent_real_roots_MeV":roots,"root_count":len(roots),
        "deltaP_low_endpoint":first,"deltaP_high_endpoint":last,
        "slope_deltaP_at_transition_fm3":slope,
        "root_difference_from_native_MeV":error,"pressure_difference_at_native_root":residual,
        "pressure_relative_residual":residual/max(abs(pt),1e-300),
        "claim":"single forward switch with no hidden recrossing on entire common support of the SAME two connected components",
        "excluded_geometry_components":"not interpolated; no global microphysical ground-state claim"}


# A cubic with three interior roots tests concealed recross detection.
assert np.allclose(unit_roots(Polynomial.fromroots([.1,.5,.9])),[.1,.5,.9],rtol=0,atol=1e-12)
manifest=load(BASE/"prepared/MANIFEST.json")
records={r["id"]:r for r in manifest["source_records"]}
campaign=load(BASE/"campaign/CAMPAIGN.json")
results=[]
for cell in campaign["cells"]:
    for i,r in enumerate(cell["candidate_results"]):
        if r["status"]!="GENERATED_VERIFIED_FINITE_DOMAIN":continue
        c,k=records[cell["crust_source_id"]],records[cell["core_source_id"]]
        lc,hc=r["crust_component"],r["core_component"]
        low=_segment_eosdata(read_prepared(c),(lc["phase_id"],lc["start_zero_based"],lc["stop_exclusive"]))
        high=_segment_eosdata(read_prepared(k),(hc["phase_id"],hc["start_zero_based"],hc["stop_exclusive"]))
        result=audit(low,high,r["transition"])
        result.update(cell=cell["id"],category=cell["category"],candidate_index=i,
                      crust_source_id=c["id"],core_source_id=k["id"])
        results.append(result)
record={"schema":"ns_cowling.new_model_independent_envelope_audit.v1",
    "algorithm":"independent SciPy CubicHermiteSpline with physical slopes n=dP/dmu; explicit normalized cubic coefficients on every union interval; all complex polynomial roots inspected",
    "script_sha256":sha(__file__),"manifest_sha256":sha(BASE/"prepared/MANIFEST.json"),
    "hidden_three_root_control":"PASS","results":results,
    "status":"PASS" if results and all(r["status"]=="PASS" for r in results) else "FAIL"}
save(BASE/"matching_validation/ENVELOPE_AUDIT.json",record)
print(json.dumps({"status":record["status"],"outputs_checked":len(results),
    "max_root_difference_MeV":max(x["root_difference_from_native_MeV"] or 0 for x in results)}))
if record["status"]!="PASS":raise SystemExit(1)
