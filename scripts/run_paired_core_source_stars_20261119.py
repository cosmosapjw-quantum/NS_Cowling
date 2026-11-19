#!/usr/bin/env python3
"""Source-original vs newly computed same-named RMF core, *same crust closure*.

Keep EXACT same frozen surface, outer and lower connector; reconstruct inner
crust with first law without crossing its original phase boundary; replace
ONLY the high-density core through the explicitly labelled native CSS or
Maxwell operation. Compare stellar observations at fixed mass, not at n_c.
Never infer that generated-core Maxwell is an original-provider coexistence.
"""
from __future__ import annotations
import hashlib,json,sys,traceback
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"scripts"),str(ROOT/"python_frontend")]
from compose_eos import read_compose,match_crust_core
from compose_eos.api import construct_maxwell_transition,_thermo_only_view
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.source_calibration import first_law_segment_audit
from compose_eos.maxwell_v2 import _segment_eosdata
from run_official_original_unified_controls_20261119 import _view
from run_fixed_mass_maxwell_css_20261119 import run_one as run_at_1p4

DATA=ROOT/"data/cowling_maxwell_native_20261119"
OLD=ROOT/"data/cowling_eos_source_calibration_20261119"
OUT=ROOT/"data/cowling_rmf_reconciliation_20261119/paired_core"
PROVIDER=DATA/"provider_original_unified"
PRIOR=json.loads((OLD/"phase_local_fixed_mass/RESULT.json").read_text())
SOURCE_SHA=json.loads((OLD/"source_audit/RESULT.json").read_text())

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def error(ex):return dict(type=type(ex).__name__,message=str(ex)[:420])

def original_data(name):
    target=PROVIDER/name/("validated_monotone_prefix_projection" if name=="DDME2" else "raw")
    original=read_compose(target)
    d=original
    ph=np.asarray(d.Iphase)
    i,j=np.flatnonzero(ph==7)[0],np.flatnonzero(ph==6)[0]
    end=int(np.searchsorted(d.nb,.49,side="right"))
    low=_view(d,int(i),int(j))
    high=_view(d,int(j),end)
    audit=first_law_segment_audit(low.nb,low.p,low.eps,low.mub)
    assert audit["valid_pressure_primary"]
    low.Q3[:]=audit["corrected_mu"]/low.mn-1
    low.Q6[:]=audit["corrected_f_energy_per_baryon"]/low.mn-1
    low.Q7[:]=low.Q6.copy()
    return low,high

def strict_baseline(name):
    p=OLD/"phase_local_repairs/whole"/name/"phasewise_firstlaw_integrated/css/REPRESENTATION.json"
    return p,json.loads(p.read_text())

def attempt(name,method):
    lower,source_high=original_data(name)
    base_path,contract=strict_baseline(name)
    if name not in ("FSU2","TM1e","DDME2"):
        raise ValueError("unrecognized same-named source model")
    high_path=DATA/"same_rmf/generated_core"/("GPPVA_"+name)
    generated=read_compose(high_path)
    # The generated physical source is NOT the provider's core.
    high=_thermo_only_view(generated)
    if method=="css":
        merged,rep=match_crust_core(lower,high,name+"_OFFICIAL_INNER",
                                   name+"_NEW_RMF_CORE","css",40,drop_composition=True)
        j=int(rep["css_transition_low_index"])
        pt=float(rep["css_transition_pressure"])
    elif method=="maxwell":
        merged,rep=construct_maxwell_transition(
             lower,high,assume_local_charge_neutrality=True,
             drop_composition=True,return_report=True)
        j=int(rep["low_phase_rows"])-1
        pt=float(rep["pressure_transition"])
    else:
        raise ValueError("unknown interface physics")
    if not rep["consistency"]["overall"]["pass_"]:
        raise RuntimeError("native finite consistency gate failed")
    n=np.asarray(merged.nb);p=np.asarray(merged.p);e=np.asarray(merged.eps)
    if not (0<=j<len(n)-1 and n[j+1]>n[j]>0 and
            np.isclose(p[j],p[j+1],atol=1e-10,rtol=1e-11)):
        raise RuntimeError("native interface physical pressure/density gate fails")
    # Stronger than keeping only outer layers: all actual inner-crust
    # original phase-7 n, P and epsilon samples to the native CSS endpoint
    # are exactly identical in the official-core and generated-core controls.
    if method=="css":
        origin=contract["upper"]
        src_j=int(origin["sharp_interfaces"][0])
        for field,values in (("n",n),("p",p),("epsilon",e)):
            a=np.asarray(origin[field][:src_j+1],dtype=float)
            b=np.asarray(values[:j+1],dtype=float)
            if not np.array_equal(a,b):
                raise RuntimeError("original inner-crust "+field+
                                   " knots changed when replacing only the core")
        original_inner_samples_unchanged=src_j+1
    else:
        original_inner_samples_unchanged=None
    lower_fingerprint={k:hashlib.sha256(json.dumps(contract[k],sort_keys=True).encode()).hexdigest()
                      for k in ("outer","lower_bridge","surface")}
    adapted=dict(surface=contract["surface"],outer=contract["outer"],
                 lower_bridge=contract["lower_bridge"],
                 upper=dict(n=list(map(float,n)),p=list(map(float,p)),
                            epsilon=list(map(float,e)),sharp_interfaces=[j]))
    whole=WholeCrustBarotrope.from_contract(adapted)
    # Both compared stars retain the SAME source crust up to transition.
    for key,fp in lower_fingerprint.items():
        got=hashlib.sha256(json.dumps(adapted[key],sort_keys=True).encode()).hexdigest()
        if got!=fp:raise RuntimeError("crust or lower bridge changed in paired control")
    if method=="css":
        ref=WholeCrustBarotrope.from_contract(contract)
        for pressure in np.geomspace(max(whole.pressure_domain[0]+1e-10,
                     whole.outer.p[0]*1.01),min(pt*.97,float(whole.outer.p[-1]*1.08)),8):
            a=whole.state_at_pressure(float(pressure))
            b=ref.state_at_pressure(float(pressure))
            if any(abs(a[q]-b[q])>1e-9*max(1.,abs(b[q])) for q in ("p","epsilon")):
                raise RuntimeError("paired low-pressure source states changed")
    dest=OUT/"contracts"/name/method/"REPRESENTATION.json"
    dump(dest,adapted)
    star=run_at_1p4(name+"_generated_"+method,dest)
    return dict(status="PAIRED_THREE_LAYER_RUN" if star.get("status","").startswith(
                 "FIXED_1P4_") else "STELLAR_SCOPE_BLOCKED",
           finite_method="PASS",method=method,source_generated=True,
           crust_frozen_against_provider=True,
           exactly_identical_original_inner_samples=original_inner_samples_unchanged,
           transition_p=pt,transition_n_low=float(n[j]),
           transition_n_high=float(n[j+1]),transition_delta_epsilon=float(e[j+1]-e[j]),
           generated_core_thermo_sha256=sha(high_path/"eos.thermo"),
           original_provider_core_sha256=sha(PROVIDER/name/"raw/eos.thermo"),
           compared_provider_core_physics="original unified Iphase=6 first-law corrected",
           existing_provider_css_1p4=PRIOR["cases"][name]["status"],
           source_density_core_difference=SOURCE_SHA["models"][name][
                "original_homogeneous_core_same_density_comparison"]["pressure_abs_rel_max"],
           source_contract_hash=sha(dest),source_contract_path=str(dest.relative_to(ROOT)),
           star=star)

def main():
    out=dict(schema="ns_cowling.stage05.paired_same_crust_source_core_reconciliation.v1",
        parent_source_audit_sha256=sha(OLD/"source_audit/RESULT.json"),
        parent_phase_local_sha256=sha(OLD/"phase_local_repairs/RESULT.json"),
        preserved_lower_layers=("outer","surface","lower_bridge"),
        claim="only upper generated RMF core and explicit matching method change; no fitted EOS",
        cases={},status="RUNNING")
    for model in ("FSU2","TM1e","DDME2"):
        entry={}
        for method in ("css","maxwell"):
            try:
                x=attempt(model,method)
            except Exception as ex:
                x=dict(status="NATIVE_OR_STELLAR_REJECTED",failure=error(ex))
            entry[method]=x
            print("CORE_PAIRED",model,method,x["status"],
                  x.get("transition_p"),x.get("star",{}).get("native_TOV",{}).get("R_km"),
                  x.get("failure"),flush=True)
        out["cases"][model]=entry
        dump(OUT/"PARTIAL.json",out)
    out["status"]="ALL_BOUNDED_PAIRED_CORE_CHOICES_TESTED"
    dump(OUT/"RESULT.json",out)
    print("PAIRED_FINISHED",list(out["cases"]),flush=True)

if __name__=="__main__":main()
