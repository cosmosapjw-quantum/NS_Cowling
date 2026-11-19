#!/usr/bin/env python3
"""Conservative source repairs at EXACT GPPVA unified inner/core boundaries.

Compare original Q3, metadata-only Euler Q3, and pressure-primary first-law
on each original phase separately. Never integrate the first law through
the *unverified* Iphase7->Iphase6 phase gap to claim a microscopic EOS.
"""
from __future__ import annotations
import hashlib,json,sys,traceback
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos import read_compose,match_crust_core
from compose_eos.api import construct_maxwell_transition
from compose_eos.maxwell_v2 import _segment_eosdata
from compose_eos.whole_crust import WholeCrustBarotrope
from compose_eos.source_calibration import first_law_segment_audit
from run_official_original_unified_controls_20261119 import _view
from build_original_fsu2_same_rmf_three_layer_20261119 import mode_result

BASE=ROOT/"data/cowling_maxwell_native_20261119"
SOURCE=BASE/"provider_original_unified"
OUT=ROOT/"data/cowling_eos_source_calibration_20261119/phase_local_repairs"
N_C=.30
SOURCE_STATUS=json.loads((SOURCE/"ORIGINAL_CONTROLS.json").read_text())

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def gate(e):return dict(type=type(e).__name__,error=str(e)[:420])

def old_source(name):
    row=SOURCE_STATUS["models"][name]
    path=(ROOT/row["source_projection"]["projection_path"]
           if row.get("source_projection") else SOURCE/name/"raw")
    d=read_compose(path)
    ph=np.asarray(d.Iphase)
    inner=np.flatnonzero(ph==7);core=np.flatnonzero(ph==6)
    if not(len(inner)>5 and len(core)>5 and inner[-1]+1==core[0]):
        raise ValueError("not exactly one original inner-core phase interface")
    i=int(inner[0]);j=int(core[0])
    end=int(np.searchsorted(d.nb,.49,side="right"))
    return d,i,j,end

def variant(low,high,law):
    l=low.copy();h=high.copy()
    before=[(np.asarray(d.nb).copy(),np.asarray(d.p).copy(),np.asarray(d.eps).copy())
            for d in (l,h)]
    if law=="raw_q3_original":
        pass
    elif law=="euler_q3_only":
        for d in (l,h):
            d.Q3[:]=(d.eps+d.p)/(d.nb*d.mn)-1.
    elif law=="phasewise_firstlaw_integrated":
        for d in (l,h):
            audit=first_law_segment_audit(d.nb,d.p,d.eps,d.mub)
            if not audit["valid_pressure_primary"]:
                raise ValueError("phasewise pressure-primary interval fails causal gate")
            d.Q3[:]=audit["corrected_mu"]/d.mn-1
            d.Q6[:]=audit["corrected_f_energy_per_baryon"]/d.mn-1
            d.Q7[:]=d.Q6.copy()
    else:
        raise ValueError("invalid source-repair policy")
    for d,(n,p,e) in zip((l,h),before):
        if not(np.array_equal(d.nb,n) and np.array_equal(d.p,p)):
            raise AssertionError("source density or pressure modified by nonregridding repair")
        if law in ("raw_q3_original","euler_q3_only") and not np.array_equal(d.eps,e):
            raise AssertionError("Q3-only metadata repair altered physical epsilon")
    return l,h

def try_method(low,high,law,model,method,outer):
    try:
        if method=="maxwell":
            eos,rep=construct_maxwell_transition(low,high,
                 assume_local_charge_neutrality=True,drop_composition=True,
                 return_report=True)
            j=int(rep["low_phase_rows"])-1
            pt=float(rep["pressure_transition"])
        else:
            eos,rep=match_crust_core(low,high,model+" INNER",model+" CORE",
                 "css",40,drop_composition=True)
            j=int(rep["css_transition_low_index"])
            pt=float(rep["css_transition_pressure"])
        if not rep["consistency"]["overall"]["pass_"]:
            raise RuntimeError("native method reports invalid overall thermodynamics")
        n,p,epsilon=np.asarray(eos.nb),np.asarray(eos.p),np.asarray(eos.eps)
        if not (0<=j<len(n)-1 and n[j+1]>n[j] and
                abs(p[j]-p[j+1])<2e-11*max(1,abs(pt))):
            raise ValueError("source transition has no compatible sharp pressure/density jump")
        match=dict(status="PASS_NATIVE_FINITE",
                  transition_pressure_MeV_fm3=pt,
                  n_low=float(n[j]),n_high=float(n[j+1]),
                  energy_jump_MeV_fm3=float(epsilon[j+1]-epsilon[j]),
                  native_consistency=True)
        # A 3-layer stellar model requires every same-phase Hermite interval
        # causal. A native finite doublet ALONE does not approve it.
        try:
            contract=dict(surface=outer["surface"],outer=outer["outer"],
              lower_bridge=outer["lower_bridge"],
              upper=dict(n=n.tolist(),p=p.tolist(),epsilon=epsilon.tolist(),
                         sharp_interfaces=[j]))
            bar=WholeCrustBarotrope.from_contract(contract)
            dest=OUT/"whole"/model/law/method/"REPRESENTATION.json"
            save(dest,contract)
            star=mode_result(bar,model,law+"_"+method,N_C)
            match["stellar"]=dict(status="NATIVE_TOV_AND_COWLING_ATTEMPTED",
                     evidence=star,contract_path=str(dest.relative_to(ROOT)))
        except Exception as ex:
            match["stellar"]=dict(status="THREE_LAYER_BLOCKED",
                                 error=gate(ex))
        return match
    except Exception as ex:
        return dict(status="NATIVE_REJECTED",error=gate(ex))

def run_model(name):
    d,i,j,end=old_source(name)
    low=_view(d,i,j)
    high=_view(d,j,end)
    reference=json.loads((SOURCE/name/"thermodynamically_corrected_smooth/REPRESENTATION.json").read_text())
    out=dict(model=name,source_hash=sha(SOURCE/name/"raw/eos.thermo"),
             physical_tiers={
               "source_original":"P,n,epsilon,Q3 unmodified; only explicit thermo-only projection",
               "euler_q3_only":"P,n,epsilon original, Q3 replaced by algebraic (epsilon+P)/n",
               "phasewise_firstlaw_integrated":"P,n original, epsilon/Q3 re-integrated *within each original phase* anchored separately, source phase seam unbridged"},
             original_phase_layout=dict(Iphase7=[i,j],Iphase6=[j,end]),
             original_mu_gap_MeV=float(high.mub[0]-low.mub[-1]),
             variants={})
    for kind in ("raw_q3_original","euler_q3_only","phasewise_firstlaw_integrated"):
        try:
            a,b=variant(low,high,kind)
            rec=dict(status="READ_ONLY_METHOD_ATTEMPT",
               original_density_pressure_preserved=True,
               raw_epsilon_preserved=kind!="phasewise_firstlaw_integrated",
               mu_gap_MeV=float(b.mub[0]-a.mub[-1]),
               firstlaw_changes_low_maxabs_MeV=float(np.max(abs(a.eps-low.eps)/a.nb)),
               firstlaw_changes_high_maxabs_MeV=float(np.max(abs(b.eps-high.eps)/b.nb)))
            for method in ("maxwell","css"):
                rec[method]=try_method(a,b,kind,name,method,reference)
            out["variants"][kind]=rec
            print("PHASE_LOCAL",name,kind,"mu_gap",rec["mu_gap_MeV"],
                "maxwell",rec["maxwell"]["status"],"css",rec["css"]["status"],
                "stellar",rec["css"].get("stellar",{}).get("status"),flush=True)
        except Exception as ex:
            out["variants"][kind]=dict(status="CASE_BLOCKED",error=gate(ex))
    return out

def main():
    out=dict(schema="ns_cowling.original_unified_phase_local_repair.v1",
      reference_source_sha256=sha(SOURCE/"ORIGINAL_CONTROLS.json"),
      success_means="finite matching and actual same-phase causal stellar acceptance separately",
      physical_phase_prohibition="no first-law smoothing through original Iphase7->6 seam",
      cases={},status="RUNNING")
    for name in ("FSU2","TM1e","DDME2"):
        try:out["cases"][name]=run_model(name)
        except Exception as ex:
            out["cases"][name]=dict(status="BLOCKED",error=gate(ex),
                traceback=traceback.format_exc()[-1100:])
        save(OUT/"PARTIAL.json",out)
    out["status"]="COMPLETE_BOUNDED_VARIANT_SCREEN"
    save(OUT/"RESULT.json",out)
    print("PHASE_LOCAL_DONE",list(out["cases"]),flush=True)
if __name__=="__main__":main()
