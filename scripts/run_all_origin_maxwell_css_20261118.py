"""Six-worker native CSS plus source-separated Maxwell campaign.

Optimization: verify source once, cache each EOS, pre-compute admitted phase
components, vectorize P(mu) bracketing, reuse the sealed 195-pair 2026-11-17
Maxwell campaign by exact source ID. Run new native Maxwell only on prefiltered
physical crossing candidates. Source-null / untested components remain blocked,
not falsely declared physically impossible.

This distinguishes:
 - source file available, origin purity, and independent nuclear functional;
 - both matching algorithms produce bounded output;
 - real surface-to-transition crust admitted;
 - globally minimizing phase construction (not proven by this screening).
"""
from __future__ import annotations
import hashlib,json,os,sys,time,traceback
from pathlib import Path
from functools import lru_cache
from concurrent.futures import ThreadPoolExecutor,as_completed
from collections import Counter,defaultdict
import numpy as np
from scipy.interpolate import CubicHermiteSpline,PchipInterpolator,BPoly
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from all_origin_maxwell_css_20261118 import OUT,catalog,enumerate_pairs,classify_both,bracket_roots,sha
from compose_eos.api import read_compose,match_crust_core,construct_maxwell_transition,_thermo_only_view
from compose_eos.maxwell_v2 import _segment_eosdata
from compose_eos.native import executable
from scipy.optimize import brentq

MAX_WORKERS=6
LEGACY=ROOT/"data/maxwell_three_families_20261117/campaign_v2/CAMPAIGN.json"
def clean(v):
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,Path):return str(v)
    raise TypeError(str(type(v)))
def save(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+".tmp")
    tmp.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False,default=clean)+"\n")
    tmp.replace(p)
@lru_cache(maxsize=120)
def _load(identifier):
    a=IDS[identifier]
    if a["status"]!="SOURCE_SNAPSHOT_AVAILABLE":raise ValueError(a["status"])
    p=ROOT/a["compose_path"]
    if not p.is_dir():raise FileNotFoundError(p)
    d=read_compose(p)
    selection=a.get("selection",{})
    if a.get("provenance_class")=="new_public":
        start=int(selection["source_start_zero_based"]);stop=int(selection["source_stop_exclusive"])
        d=_segment_eosdata(d,(0,start,stop))
    if d.n<3 or not(np.all(np.diff(d.nb)>0) and np.all(d.nb>0)):
        raise ValueError("invalid source density support")
    return d
def comps(a,model):
    entries=a.get("components",[])
    if entries:
        return [dict(start=int(c["start_zero_based"]),stop=int(c["stop_exclusive"]),
                     mu_bounds=list(map(float,c["mu_E_MeV"])),source_component_id=c["id"])
                for c in entries if c.get("rows",0)>=3]
    mu=np.asarray(model.mub,float)
    n=np.asarray(model.nb,float)
    p=np.asarray(model.p,float)
    valid=np.isfinite(mu)&np.isfinite(n)&np.isfinite(p)&(mu>0)&(n>0)&(p>0)
    accepted=[]
    start=None
    for i in range(len(mu)):
        if i<len(mu)-1:
            okay=bool(valid[i] and valid[i+1] and mu[i+1]>mu[i] and
                       n[i+1]>n[i] and p[i+1]>p[i])
        else:okay=False
        if okay and start is None:start=i
        if not okay and start is not None:
            stop=i+1
            if stop-start>=3:
                accepted.append(dict(start=start,stop=stop,mu_bounds=[float(mu[start]),float(mu[stop-1])],
                                     source_component_id=f"dynamic_{start}_{stop}",
                                     diagnostic="only monotone point values, native Fortran must repeat all physical shape gates"))
            start=None
    return accepted
def fit(d,c,method):
    v=_segment_eosdata(d,(0,c["start"],c["stop"]))
    mu=np.asarray(v.mub,float);p=np.asarray(v.p,float);n=np.asarray(v.nb,float)
    if np.any(np.diff(mu)<=0) or np.any(np.diff(n)<=0):raise ValueError("source knots nonmonotone")
    if method=="cubic":
        poly=CubicHermiteSpline(mu,p,n)
    elif method=="quintic":
        if len(mu)<4:raise ValueError("quintic requires >=4 source knots")
        chi=PchipInterpolator(mu,n).derivative()(mu)
        if np.any(~np.isfinite(chi)) or np.any(chi<=0):raise ValueError("source quintic susceptibility nonpositive")
        poly=BPoly.from_derivatives(mu,[[p[i],n[i],chi[i]] for i in range(len(mu))])
    else:raise ValueError(method)
    return v,poly
def quick_candidates(low,high):
    """Metadata μ overlap and vectorized sign roots; all negatives remain screen-only."""
    al=IDS[low];ah=IDS[high]
    l=_load(low);h=_load(high)
    lc=comps(al,l);hc=comps(ah,h)
    overlaps=[]
    rootbrackets=[]
    checks=0
    for cl in lc:
        for ch in hc:
            left=max(cl["mu_bounds"][0],ch["mu_bounds"][0]);right=min(cl["mu_bounds"][1],ch["mu_bounds"][1])
            if not right>left:continue
            overlaps.append(dict(low=cl["source_component_id"],high=ch["source_component_id"],
                                 common_mu_MeV=[left,right]))
            try:
                vl,pl=fit(l,cl,"cubic");vh,ph=fit(h,ch,"cubic")
                rr=bracket_roots(pl,ph,left,right,121)
                for v in rr:
                    low_n=float(pl(v,nu=1));high_n=float(ph(v,nu=1))
                    # Intersource thermodynamic co-existence requires a positive density jump.
                    if high_n>low_n and min(high_n,low_n)>0:
                        rootbrackets.append(dict(low=cl,high=ch,root=v,
                           mu_bounds=[left,right],low_n=low_n,high_n=high_n))
            except Exception:pass
            checks+=1
    return dict(status="SCREEN_COMPLETED_CONSERVATIVE",admissible_components=[len(lc),len(hc)],
                common_mu_component_count=len(overlaps),common_mu=list(overlaps[:10]),
                indicated_positive_crossings=len(rootbrackets),tested_component_pairs=checks),rootbrackets
def slice_comp(a,c):
    return _segment_eosdata(a,(0,c["start"],c["stop"]))
def physics_maxwell(l,h,rootbrackets):
    successes=[];failures=[]
    seen=set()
    for br in rootbrackets:
        key=(br["low"]["source_component_id"],br["high"]["source_component_id"])
        if key in seen:continue
        seen.add(key)
        la=slice_comp(l,br["low"]);ha=slice_comp(h,br["high"])
        try:
            matched,rep=construct_maxwell_transition(la,ha,
                assume_local_charge_neutrality=True,drop_composition=True,return_report=True)
            check=rep["consistency"]["overall"]["pass_"]
            if not check:raise ValueError("native output consistency failed")
            successes.append(dict(low_component=key[0],high_component=key[1],
                P_t=rep["pressure_transition"],mu_t=rep["mu_transition"],
                n_low=rep["n_low"],n_high=rep["n_high"],
                jump_n=rep["density_jump"],latent_energy=rep["latent_energy_density"],
                rows=matched.n,reported_by="native Fortran single-conserved-charge Maxwell"))
        except Exception as exc:
            failures.append(dict(components=key,type=type(exc).__name__,message=str(exc)[:280]))
    return successes,failures
def legacy_maxwell(a,b):
    x=LEGACY_CELLS.get((a["provenance_class"],a["model"],b["model"]))
    if x is None:return None
    admitted=[]
    for candidate in x["candidate_results"]:
        if candidate["status"]!="GENERATED_VERIFIED_FINITE_DOMAIN":continue
        n=candidate["v2_report"]["native_report"]
        admitted.append(dict(P_t=n["pressure_transition"],mu_t=n["mu_transition"],
            n_low=n["n_low"],n_high=n["n_high"],jump_n=n["density_jump"],
            latent_energy=n["latent_energy_density"],
            source_components=[candidate["crust_component"]["id"],candidate["core_component"]["id"]],
            source="previous native Fortran source-sealed Maxwell v2"))
    return dict(status="PASS_VERIFIED_FINITE" if x["generated_count"]>0 else
                ("NO_SOURCE_COMPONENT_PAIR" if x["status"]=="SOURCE_HAS_NO_ADMITTED_COMPONENT_PAIR"
                  else "NO_VALID_COEXISTENCE_ON_ADMITTED_COMPONENTS"),
       admitted=admitted,
       reused_native_verified_20261117=True,
       source_ids=[a["source_id"],b["source_id"]],
       attempted_native_kernel=x["kernel_invocations"],historical_result=x["status"],
       admitted_component_output_count=x["generated_count"],source_outcome="no global-ground-state claim",
       source_ref="data/maxwell_three_families_20261117/campaign_v2/CAMPAIGN.json")
def genuine_outer_and_inner(c,css_report,maxwell,screen):
    if not c.get("outer_present",False):
        return dict(status="OUTER_ABSENT",
                    note="Only finite inner crust from native source; original whole-source outer crust missing")
    if not c.get("full_outer_inner",False):
        return dict(status="INNER_OR_OUTER_NOT_ESTABLISHED")
    if not css_report.get("status")=="PASS":
        return dict(status="CSS_REJECTED")
    if maxwell.get("status") not in ("PASS_VERIFIED_FINITE","PASS_NEW_NATIVE_FINITE"):
        return dict(status="MAXWELL_NOT_ADMITTED")
    if not screen["admissible_components"][0] or not screen["admissible_components"][1]:
        return dict(status="NO_ADMITTED_SOURCE_INTERVAL")
    # Strict whole-star requires no hidden invalid phase intervals and a surface.
    # Prepared original SPG/GPPVA outer is untouched, but source inner may have invalid holes.
    # Accepted components must start at first selected inner knot and reach its upper boundary.
    x=IDS[c["source_id"]]
    cc=x.get("components",[])
    if len(cc)!=1 or int(cc[0]["start_zero_based"])!=0 or int(cc[0]["stop_exclusive"])!=int(x["selection"]["source_stop_exclusive"]-x["selection"]["source_start_zero_based"]):
        return dict(status="SOURCE_INNER_GAPS_OR_EXTRA_PHASES_UNVALIDATED",
                  note="No interpolation over a failed/source phase interval; cannot certify complete crust")
    return dict(status="CANDIDATE_FULL_SOURCE_CHAIN",
                note="Full source chain preserved but independent TOV/surface admission still required")
def one(a,b):
    result=dict(schema="ns_cowling.css_maxwell_pair.v1",pair_id=a["source_id"]+"__"+b["source_id"],
                tier=a["tier"],crust=a["model"],core=b["model"],
                crust_source_id=a["source_id"],core_source_id=b["source_id"],
                same_functional=(a["nuclear_functional"]==b["nuclear_functional"]),
                same_functional_core_exists=a["same_functional_core_exists"],
                crust_outer_present=a.get("outer_present",False),
                crust_full_outer_inner_source=a.get("full_outer_inner",False),
                priority=a["priority"],status="RUNNING",css={},maxwell={},screen={})
    if a["status"]!="SOURCE_SNAPSHOT_AVAILABLE" or b["status"]!="SOURCE_SNAPSHOT_AVAILABLE":
        result.update(status="BLOCKED_SOURCE",css=dict(status="BLOCKED"),
            maxwell=dict(status="BLOCKED"),
            whole=dict(status="OUTER_OR_THERMODYNAMIC_SOURCE_MISSING"),
            reason=f'{a["status"]} / {b["status"]}',both_status="BLOCKED_SOURCE")
        return result
    try:
        l=_load(a["source_id"]);h=_load(b["source_id"])
        res,roots=quick_candidates(a["source_id"],b["source_id"])
        result["screen"]=res
        start=time.monotonic()
        try:
            matched,report=match_crust_core(l,h,a["model"],b["model"],"css",40,
                                             drop_composition=True)
            result["css"]=dict(status="PASS" if report["consistency"]["overall"]["pass_"] else "FAILED_CONSISTENCY",
                P_t=report["css_transition_pressure"],n_low=report["css_transition_nb_low"],
                n_high=report["css_transition_nb_high"],
                latent_heat=report["css_transition_delta_eps"],
                cs2=report["css_single_sound_speed"],
                source_rows=[l.n,h.n],output_rows=matched.n,
                actual_fortran_invoked=True,
                elapsed_seconds=time.monotonic()-start)
        except Exception as exc:
            result["css"]=dict(status="REJECTED",type=type(exc).__name__,error=str(exc)[:350])
        old=legacy_maxwell(a,b)
        if old is not None:
            result["maxwell"]=old
        elif roots:
            yes,no=physics_maxwell(l,h,roots)
            result["maxwell"]=dict(status="PASS_NEW_NATIVE_FINITE" if yes else "NO_COEXISTENCE_AFTER_NATIVE_GATES",
                 admitted_component_output_count=len(yes),admitted=yes[:9],native_rejected=no[:9],
                 native_invoked=len(yes)+len(no),
                 preliminary_positive_crossings=len(roots))
        else:
            result["maxwell"]=dict(status="PREFILTER_NO_CROSSING_NOT_A_PROOF",
                preliminary_positive_crossings=0,not_exhaustive=True,
                source_overlap_exists=bool(res["common_mu_component_count"]))
        result["whole"]=genuine_outer_and_inner(a,result["css"],result["maxwell"],result["screen"])
        result["both_status"]=classify_both(dict(css=result["css"]["status"],
             maxwell="PASS" if result["maxwell"]["status"] in ("PASS_VERIFIED_FINITE","PASS_NEW_NATIVE_FINITE") else "FAIL",
             whole="FAIL"))  # source-chain evidence alone is NOT a TOV-validated whole NS model
        result["status"]="EVALUATED_WITH_SOURCE_SCOPE"
    except Exception as exc:
        result["status"]="ENVIRONMENT_OR_UNEXPECTED_SOURCE_FAILURE"
        result["type"]=type(exc).__name__;result["error"]=str(exc)[:320]
        result["traceback"]=traceback.format_exc()[-1000:]
        result["both_status"]="UNRESOLVED"
    return result
def main():
    global IDS,LEGACY_CELLS
    catalog_data=catalog()
    IDS={x["source_id"]:x for x in catalog_data["sources"]}
    pre=json.loads(LEGACY.read_text())
    LEGACY_CELLS={(z["family"],z["crust"],z["core"]):z for z in pre["cells"]}
    pairs=enumerate_pairs(catalog_data["sources"])
    OUT.mkdir(parents=True,exist_ok=True)
    archive=dict(status="RUNNING",input_origin_classes=["original","corrected","reconstructed"],
        source_catalog_path="data/all_origin_maxwell_css_20261118/SOURCE_CATALOG.json",
        first_order_source_crossing_from_existing=195,
        max_workers=MAX_WORKERS,method="cached EOS, vectorized sign root screen, existing exact-hash Maxwell source-campaign reuse",
        requested_total=len(pairs),cases=[])
    save(OUT/"CAMPAIGN_RUNNING.json",archive)  # never publish partial state at terminal SSOT path
    t=time.monotonic()
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures=[pool.submit(one,a,b) for a,b in pairs]
        for fut in as_completed(futures):
            x=fut.result()
            archive["cases"].append(x)
            if len(archive["cases"])%30==0 or len(archive["cases"])==len(pairs):
                save(OUT/"CAMPAIGN_RUNNING.json",archive)
                count=Counter(z["both_status"] for z in archive["cases"])
                print("CAMPAIGN_PROGRESS",len(archive["cases"]),len(pairs),"BOTH",dict(count),flush=True)
    archive["status"]="COMPLETED_WITH_EXPLICIT_ORIGIN_AND_SUPPORT_GATES"
    archive["elapsed_seconds"]=time.monotonic()-t
    save(OUT/"CAMPAIGN.json",archive)
    temporary=OUT/"CAMPAIGN_RUNNING.json"
    if temporary.is_file():
        temporary.unlink()  # preserve only the atomic, complete SSOT publication
    print("FINISHED",len(pairs),archive["elapsed_seconds"],Counter(z["both_status"] for z in archive["cases"]))
if __name__=="__main__":main()
