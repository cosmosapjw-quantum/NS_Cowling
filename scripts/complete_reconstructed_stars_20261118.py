"""Link strictly origin-separated reconstructed crust/core candidates to the
reconstructed BSk22 outer plus anchored surface of the prior E worktree.

Result is phenomenological COMPLETE STAR, not a globally unified microphysical
crust; every fitted lower bridge and chosen join index are exposed. Never
reclassify an unverified original-CompOSE source or external nuda input.
"""
from __future__ import annotations
from pathlib import Path
import hashlib,json,sys,time,traceback
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from compose_eos import read_compose,match_crust_core
from compose_eos.api import construct_maxwell_transition
from compose_eos.whole_crust import WholeCrustBarotrope,CausalExponentialInterval
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS,solve_sharp_background
from all_origin_maxwell_css_20261118 import OUT,catalog,sha
SOURCE_E=ROOT/"data/stellar_extension_20261117/eos/fine/REPRESENTATION.json"
INPUT=OUT/"CAMPAIGN.json"
OUTPUT=OUT/"whole_stars"
CASES=("SPG_M1","SPG_M2","SPG_M3","SPG_M4","SPG_M5","GPPVA_FSU2")
def save(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def eos_source(name,role):
    if role=="crust":
        if name.startswith("SPG_"):
            return ROOT/f"data/research_matching_20261117/spg/native_129_repaired/{name}/generated"
        if name=="GPPVA_FSU2":
            return ROOT/"data/research_native_20261116/source_grid_cohort/GPPVA_FSU2/attempt_001/output/generated"
        raise ValueError(name)
    if name=="PK1_CANONICAL":
        return ROOT/"data/new_model_maxwell_20261117/generated/PK1_CANONICAL/fine/compose"
    if name=="GM_GM1":
        return ROOT/"data/maxwell_three_families_20261117/native_core_extended_v2/GM_GM1/generated"
    raise ValueError(name)
def fixed_mass(e):
    bg=SharpContractEnthalpyEOS(e)
    controls=dict(r0_km=.001,rtol=2e-9,atol=2e-11,max_step_h=.004)
    def one(n):return solve_sharp_background(bg,n,**controls)
    left=.21;right=.53
    if not one(left).mass_msun<1.4<one(right).mass_msun:
        raise ValueError("1.4 Msun not bracketed")
    root=brentq(lambda n:one(n).mass_msun-1.4,left,right,xtol=1e-9)
    res=one(root)
    if abs(res.mass_msun-1.4)>2e-7:raise ValueError("mass root residual not admitted")
    return dict(mass_msun=res.mass_msun,R_km=res.radius_km,
      central_density_fm3=res.central_density_fm3,
      radius_interface=[x["r_km"] for x in res.interfaces],
      interface_count=len(res.interfaces))
def try_whole(c,co,method,reference):
    if method=="maxwell":
        matched,info=construct_maxwell_transition(c,co,
            assume_local_charge_neutrality=True,drop_composition=True,return_report=True)
        j=int(info["low_phase_rows"])-1
        pt=float(info["pressure_transition"])
    elif method=="css":
        matched,info=match_crust_core(c,co,"crust","core","css",40,drop_composition=True)
        j=int(info["css_transition_low_index"])
        pt=float(info["css_transition_pressure"])
    else:raise ValueError(method)
    n=np.asarray(matched.nb,float);p=np.asarray(matched.p,float);e=np.asarray(matched.eps,float)
    if not 0<=j<len(n)-1:raise ValueError("missing sharp jump label")
    if not abs(p[j]-p[j+1])<2e-11*max(1,abs(pt)):
        raise ValueError("sharp pressure mismatch")
    if not n[j+1]>n[j]:raise ValueError("wrong density jump orientation")
    surface=reference["surface"];outer=reference["outer"]
    end=tuple(map(float,(outer["n_fm3"][-1],outer["p_MeV_fm3"][-1],
        outer["n_fm3"][-1]*outer["F_MeV"][-1])))
    admitted=None
    rejections=[]
    # Cheap deterministic finite scan, avoids all arbitrary bridge shifts.
    for k in range(min(50,j-2)):
        if n[k]<=end[0]:continue
        right=tuple(float(t[k]) for t in (n,p,e))
        try:
            b=CausalExponentialInterval(*end,*right)
            status=b.validate()
            if status["status"]!="PASS":raise ValueError("connector causal gate")
            representation=dict(surface=surface,outer=outer,lower_bridge=b.contract(),
                upper=dict(n=n[k:].tolist(),p=p[k:].tolist(),epsilon=e[k:].tolist(),
                           sharp_interfaces=[j-k]))
            whole=WholeCrustBarotrope.from_contract(representation)
            star=fixed_mass(whole)
            admitted=dict(join_index=k,join_n_fm3=float(n[k]),join_pressure=float(p[k]),
                lower_bridge_kind="causal_exponential_reconstructed",
                lower_bridge_validation=status,upper_source_unchanged=True,
                stellar=star,absolute_transition_index=j,transition_P=pt,
                density_jump=float(n[j+1]-n[j]),
                latent_heat=float(e[j+1]-e[j]),
                source_n_min=float(n[0]),source_n_max=float(n[-1]),
                output_provider="existing code reconstructed E BSk22 lower/outer/surface")
            break
        except Exception as exc:
            rejections.append(dict(index=k,category=type(exc).__name__,why=str(exc)[:140]))
    return dict(status="CONDITIONAL_WHOLE_NS_VERIFIED_TOV" if admitted else "WHOLE_NS_NOT_ADMITTED",
                actual_method=method,matched_source_rows=matched.n,
                match_P_t=pt,matched_transition_j=j,admitted=admitted,
                native_report={
                    "Maxwell_pressure_transition":info.get("pressure_transition"),
                    "CSS_interface_pressure":info.get("css_transition_pressure"),
                    "source_matched_consistency":info.get("consistency",{}).get("overall",{}).get("pass_")},
                rejections=rejections[:8])
def worker(name):
    b=next(c for c in GLOBAL["cases"] if c["crust"]==name and c["tier"]=="reconstructed" and
           c["both_status"]=="FINITE_ONLY_NOT_WHOLE")
    core=b["core"]
    folder=OUTPUT/name
    folder.mkdir(parents=True,exist_ok=True)
    l=read_compose(eos_source(name,"crust"))
    h=read_compose(eos_source(core,"core"))
    ref=json.loads(SOURCE_E.read_text())
    outcome=dict(source=name,core=core,
       nuclear_functional_distinct=b["same_functional"] is False,
       outer_reference_sha256=sha(SOURCE_E),
       source_low_sha256=sha(eos_source(name,"crust")/"eos.thermo"),
       source_high_sha256=sha(eos_source(core,"core")/"eos.thermo"),
       class_scope="reconstructed code-only, no original/corrected provider mixing",
       physical_outer_scope="BSk22 outer + an anchored extrapolating surface polytrope from E, NOT source-native micro EOS",
       methods={})
    for method in ("maxwell","css"):
        try:
            outcome["methods"][method]=try_whole(l,h,method,ref)
        except Exception as exc:
            outcome["methods"][method]=dict(status="BLOCKED_OR_REJECTED",
                 type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc()[-1000:])
        save(folder/"RESULT.json",outcome)
        print("WHOLE",name,core,method,outcome["methods"][method]["status"],flush=True)
    outcome["both_conditional_whole"]=all(x["status"]=="CONDITIONAL_WHOLE_NS_VERIFIED_TOV"
                                       for x in outcome["methods"].values())
    save(folder/"RESULT.json",outcome)
    return outcome
def main():
    global GLOBAL
    GLOBAL=json.loads(INPUT.read_text())
    with ThreadPoolExecutor(max_workers=3) as pool:
        results=[z.result() for z in as_completed([pool.submit(worker,name) for name in CASES])]
    save(OUTPUT/"CAMPAIGN.json",dict(status="COMPLETE",n=len(results),
        both_conditional_whole=sum(x["both_conditional_whole"] for x in results),
        source_model_ids=[x["source"] for x in results],
        result_files=[f"{x['source']}/RESULT.json" for x in results]))
    print("WHOLE_COMPLETED",len(results),sum(x["both_conditional_whole"] for x in results))
if __name__=="__main__":main()
