#!/usr/bin/env python3
"""BBP/BPS secondary-literature numerical table + original CompOSE BPS outer.

This does not call a secondary digitization a sealed original BBP EOS.
Retain source P,epsilon,n in a separate raw branch and quantify the bounded
thermodynamic reconstruction needed for a causal three-layer model.
"""
from __future__ import annotations
import hashlib,json,re,subprocess,sys,time,traceback
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import CubicHermiteSpline,PchipInterpolator,BPoly

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos import read_compose,match_crust_core
from compose_eos.api import construct_maxwell_transition,_thermo_only_view
from compose_eos.eos_data import EOSData
from compose_eos.whole_crust import PiecewisePolytrope,SurfacePolytrope,CausalExponentialInterval,WholeCrustBarotrope
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from build_original_fsu2_same_rmf_three_layer_20261119 import mode_result

DATA=ROOT/"data/cowling_maxwell_native_20261119"
OUT=DATA/"bbp_bps"
PDF=OUT/"source_reference/diener_2008_thesis.pdf"
SRC_PDF_SHA="2f12c3c6208709ace642e30b484cb56d1b86b7de48bd6540598165a1f8f3803d"
BPS=ROOT/"data/maxwell_three_families_20261117/intake/references/APR_APR"
BBP_PHYSICAL_CRUST_MAX_N=.08
CORES={
    "GPPVA_FSU2":DATA/"same_rmf/generated_core/GPPVA_FSU2",
    "GPPVA_TM1e":DATA/"same_rmf/generated_core/GPPVA_TM1e",
    "GPPVA_DDME2":DATA/"same_rmf/generated_core/GPPVA_DDME2",
    "GM_GM1":ROOT/"data/maxwell_three_families_20261117/native_core_extended_v2/GM_GM1/generated",
    "PK1_CANONICAL":ROOT/"data/new_model_maxwell_20261117/generated/PK1_CANONICAL/fine/compose",
}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def failure(e):return dict(type=type(e).__name__,error=str(e)[:490])

def read_bbp_secondary():
    if not PDF.is_file():
        raise FileNotFoundError("Diener 2008 BBP numeric table PDF missing; download https://arxiv.org/pdf/0806.0747")
    actual=sha(PDF)
    if actual!=SRC_PDF_SHA:
        raise RuntimeError("BBP source PDF byte identity changed; re-audit Table 10.4")
    source=subprocess.run(["pdftotext","-f","109","-l","109","-layout",str(PDF),"-"],
         capture_output=True,text=True,check=True).stdout
    pattern=re.compile(r'^\s*(\d+(?:\.\d+)?)\s*×\s*10(\d{2})\s+(\d+(?:\.\d+)?)\s*×\s*10(\d{2})\s+(\d+(?:\.\d+)?)\s*×\s*10(\d{2})',re.M)
    entries=[(float(a)*10**int(b),float(c)*10**int(d),float(e)*10**int(f)) for a,b,c,d,e,f in pattern.findall(source)]
    if len(entries)!=35 or entries[0][2]!=2.67e35:
        raise ValueError("Diener Table 10.4 extraction fails source count or first row")
    # epsilon is specified as mass-equivalent g/cm3 in the secondary thesis.
    mass_g_cm3=np.array([v[0] for v in entries],float)
    dyn_cm2=np.array([v[1] for v in entries],float)
    nb_cm3=np.array([v[2] for v in entries],float)
    c=299792458.
    unit=1.602176634e32 # J/m3 per MeV/fm3
    pressure=dyn_cm2*.1/unit
    energy=mass_g_cm3*(1000*c*c)/unit
    density=nb_cm3*1e-39
    if not (np.all(np.diff(density)>0) and np.all(np.diff(pressure)>0)
            and np.all(np.diff(energy)>0)):
        raise ValueError("BBP extracted source has non-monotone raw coordinates")
    return density,pressure,energy,dict(
       source_kind="SECONDARY_LITERATURE_DIGITIZATION_FROM_BBP_REFERENCE",
       source="Diener, arXiv:0806.0747v2, Table 10.4 p97, referencing original Baym-Bethe-Pethick EOS",
       source_url="https://arxiv.org/abs/0806.0747",
       source_pdf_sha256=actual,
       original_primary_BBP_bytes="NOT_IN_HAND",
       numeric_rows=35,precision="published ~3 significant digits",
       mass_density_g_cm3_to_MeV_fm3=1000*c*c/unit,
       dyn_cm2_to_MeV_fm3=.1/unit)

def simple_eos(n,p,eps,*,mn=939.,phase=0):
    n=np.asarray(n,float);p=np.asarray(p,float);eps=np.asarray(eps,float)
    zero=np.zeros(len(n))
    mu=(eps+p)/n
    return EOSData(nb=n,mn=mn,mp=mn,Q1=p/n,Q2=zero.copy(),Q3=mu/mn-1,
        Q4=zero.copy(),Q5=zero.copy(),Q6=eps/(n*mn)-1,
        Q7=eps/(n*mn)-1,q_extra=zero.copy(),
        Iphase=np.full(len(n),phase,dtype=int),
        Yn=zero.copy(),Yp=zero.copy(),Ye=zero.copy(),Ymu=zero.copy(),
        Aav=zero.copy(),Zav=zero.copy(),Yh=zero.copy(),
        composition_file_present=False,composition_available=False,
        beta_equilibrated=True)

def make_bridge(n,p,e,bps):
    """Minimal published-rounding n(first) correction imposed by cold first law.

    Alter only the first *derived* BBP density within its 3-significant-digit
    interval. Require exact P and E as printed and a causal connector to the
    last strictly lower-P original BPS point. Nothing changes the raw table.
    """
    kk=np.flatnonzero((bps.nb<n[0]) & (bps.p<p[0]))
    if not len(kk):raise RuntimeError("no BPS original point below BBP neutron-drip pressure")
    k=int(kk[-1])
    outer=PiecewisePolytrope(bps.nb[:k+1],bps.p[:k+1],bps.eps[0]/bps.nb[0])
    assert outer.validate()["status"]=="PASS"
    last=(float(outer.nb[-1]),float(outer.p[-1]),float(outer.nb[-1]*outer.F[-1]))
    # Solve rounded first density for mu_right-mu_left ~= integral dP/n.
    mul=(last[1]+last[2])/last[0]
    a,b=n[0]-5e-7,n[0]+5e-7
    def discrepancy(nr):
        mur=(e[0]+p[0])/nr
        expected=mul+(p[0]-last[1])/(.5*(last[0]+nr))
        return mur-expected
    if discrepancy(a)*discrepancy(b)>0:
        return None,dict(status="BLOCKED_BBP_BPS_ROUNDING_INTERVAL_DOES_NOT_SUPPORT_C1_JOIN",
                         first_published_n=float(n[0]),admitted_n_range=[a,b],
                         first_mu_BBP_raw=float((e[0]+p[0])/n[0]),
                         outer_mu_end=mul,
                         residual_at_density_rounding_endpoints=[discrepancy(a),discrepancy(b)])
    n_adjust=brentq(discrepancy,a,b,xtol=1e-15)
    nr=n.copy();nr[0]=n_adjust
    try:
        bridge=CausalExponentialInterval(*last,float(n_adjust),float(p[0]),float(e[0]))
        proof=bridge.validate()
    except Exception as ex:
        return None,dict(status="BLOCKED_CAUSAL_BRIDGE",error=failure(ex))
    if proof["status"]!="PASS":
        return None,dict(status="BLOCKED_CAUSAL_BRIDGE",proof=proof)
    return (nr,outer,bridge),dict(status="ROUNDED_FIRST_N_WITHIN_PUBLISHED_PRECISION_AND_CAUSAL",
        original_first_n=float(n[0]),effective_first_n=float(n_adjust),
        delta_n_fm3=float(n_adjust-n[0]),
        max_rounding_half_width_n_fm3=5e-7,
        published_first_p_epsilon_unchanged=True,
        bps_last_source_index=k,
        bps_last_n=float(last[0]),bps_last_p=float(last[1]),
        bps_source_selected_rows=k+1,
        first_bbP_pressure=float(p[0]),
        effective_mu_BBP=float((e[0]+p[0])/n_adjust),
        source_first_mu_BBP=float((e[0]+p[0])/n[0]),
        causal_bridge_certificate=proof)

def screen_roots(low,high,order):
    m0=(np.asarray(low.eps)+np.asarray(low.p))/np.asarray(low.nb)
    m1=(np.asarray(high.eps)+np.asarray(high.p))/np.asarray(high.nb)
    def poly(d,mu,order):
        mu=np.asarray(mu,float);n=np.asarray(d.nb,float);p=np.asarray(d.p,float)
        if np.any(np.diff(mu)<=0) or np.any(np.diff(p)<=0):
            raise ValueError("nonmonotone original mu_E or pressure; cannot fit P(mu)")
        if order=="cubic":return CubicHermiteSpline(mu,p,n)
        c=PchipInterpolator(mu,n).derivative()(mu)
        if np.any(c<=0):raise ValueError("quintic source compressibility not positive")
        return BPoly.from_derivatives(mu,[[p[i],n[i],c[i]] for i in range(len(mu))])
    lo=poly(low,m0,order);hi=poly(high,m1,order)
    a=max(m0[0],m1[0]);b=min(m0[-1],m1[-1])
    if b<=a:return dict(status="NO_COMMON_MU_SUPPORT")
    xs=np.unique(np.r_[m0[(m0>=a)&(m0<=b)],m1[(m1>=a)&(m1<=b)],np.linspace(a,b,501)])
    vals=hi(xs)-lo(xs)
    found=[]
    for i in range(len(xs)-1):
        if vals[i]*vals[i+1]<0:
            mu=float(brentq(lambda t:float(hi(t)-lo(t)),xs[i],xs[i+1],xtol=1e-10))
            nl=float(lo(mu,nu=1));nh=float(hi(mu,nu=1))
            if nh>nl>0 and lo(mu)>0:
                found.append(dict(mu=mu,p=float(lo(mu)),n_low=nl,n_high=nh))
    return dict(status="SCREEN_POSITIVE_JUMP" if found else "SCREEN_NEGATIVE_NOT_NO_GO",
                roots=found,common_mu=[float(a),float(b)],test_points=len(xs))

def attempt(low,core):
    out={}
    for order in ("cubic","quintic"):
        try:out["screen_"+order]=screen_roots(low,core,order)
        except Exception as ex:out["screen_"+order]=dict(status="SCREEN_NOT_ADMISSIBLE",error=failure(ex))
    for meth in ("maxwell","css"):
        try:
            if meth=="maxwell":
                matched,rep=construct_maxwell_transition(_thermo_only_view(low),
                          _thermo_only_view(core),assume_local_charge_neutrality=True,
                          drop_composition=True,return_report=True)
                j=int(rep["low_phase_rows"])-1
                t=dict(status="PASS_NATIVE_FINITE",P_t=float(rep["pressure_transition"]),
                       mu_t=float(rep["mu_transition"]),n_low=float(rep["n_low"]),n_high=float(rep["n_high"]))
            else:
                matched,rep=match_crust_core(low,core,"BBP","CORE","css",40,
                                             drop_composition=True)
                j=int(rep["css_transition_low_index"])
                t=dict(status="PASS_NATIVE_FINITE",P_t=float(rep["css_transition_pressure"]),
                       n_low=float(rep["css_transition_nb_low"]),
                       n_high=float(rep["css_transition_nb_high"]))
            if not rep["consistency"]["overall"]["pass_"]:raise ValueError("native consistency gate failed")
            t.update(transition_index=j,rows=matched.n)
            out[meth]=t
            out[meth+"_model"]=(matched,j)
        except Exception as ex:
            out[meth]=dict(status="REJECTED",error=failure(ex))
    return out

def main():
    OUT.mkdir(exist_ok=True,parents=True)
    n,p,e,source=read_bbp_secondary()
    bps=read_compose(BPS)
    bbp_start=float(n[0])
    below=np.searchsorted(n,BBP_PHYSICAL_CRUST_MAX_N,"right")
    if below<8:raise RuntimeError("insufficient BBP inner-crust points")
    selection=dict(BBP_inner_crust_source_rows=int(below),n_max_fm3=float(n[below-1]),
         selection="fixed n<=0.08 fm^-3 independent of chosen core")
    raw=simple_eos(n[:below],p[:below],e[:below])
    matched_lower,connector=make_bridge(n[:below],p[:below],e[:below],bps)
    entries=dict(raw_BBP=raw)
    summary=dict(schema="ns_cowling.bbp_bps_original_source_and_corrected.v1",
         source=source,bps_source=dict(compose_APR_eos_path=str(BPS.relative_to(ROOT)),
         original_thermo_sha256=sha(BPS/"eos.thermo"),
         BPS_outside_unified_apr="BPS outer + SLy4 inner pre-attached; only n<=neutron drip used"),
         BBP_selection=selection,BBP_BPS_transition_contract=connector,
         source_reconstruction={},
         combinations=[],status="RUNNING")
    if matched_lower is not None:
        nr,outer,bridge=matched_lower
        bp=PiecewisePolytrope(nr,p[:below],e[0]/nr[0])
        validation=bp.validate()
        if validation["status"]!="PASS":raise RuntimeError("first-law corrected BBP fails causality")
        corrected=simple_eos(nr,p[:below],nr*bp.F)
        entries["pressure_primary_first_law_BBP"]=corrected
        summary["source_reconstruction"]=dict(
            max_abs_delta_energy_per_baryon_MeV=float(np.max(abs(bp.F-e[:below]/nr))),
            max_abs_delta_energy_density_MeV_fm3=float(np.max(abs(bp.F*nr-e[:below]))),
            original_P_unchanged=True,all_original_n_unchanged_except_first_within_rounding=True,
            original_energies_unchanged_in_raw_branch=True,
            reconstructed_thermo_status="FIRST_LAW_ANALYTIC_PIECEWISE_PRESSURE",
            max_cs2=validation["max_cs2_over_c2"])
    for branch,low in entries.items():
        for name,source_path in CORES.items():
            rec=dict(branch=branch,core=name,source_compose=str(source_path.relative_to(ROOT)),
                     source_sha256=sha(source_path/"eos.thermo"),
                     status="RUNNING")
            try:
                core=read_compose(source_path)
                result=attempt(low,core)
                rec.update(screen_cubic=result["screen_cubic"],
                     screen_quintic=result["screen_quintic"],
                     maxwell=result["maxwell"],css=result["css"])
                rec["status"]="FINITE_TEST_RECORDED"
                if matched_lower is not None and branch=="pressure_primary_first_law_BBP":
                    for method in ("maxwell","css"):
                        if result[method]["status"]!="PASS_NATIVE_FINITE":continue
                        matched,j=result[method+"_model"]
                        try:
                            surface=SurfacePolytrope.from_anchor(
                                outer.nb[0],outer.p[0],outer.nb[0]*outer.F[0])
                            contract=dict(surface=surface.contract(),
                              outer=outer.contract(),lower_bridge=bridge.contract(),
                              upper=dict(n=matched.nb.tolist(),p=matched.p.tolist(),
                                epsilon=matched.eps.tolist(),sharp_interfaces=[j]))
                            whole=WholeCrustBarotrope.from_contract(contract)
                            dest=OUT/"whole_models"/name/method
                            save(dest/"REPRESENTATION.json",contract)
                            results=mode_result(whole,"BBP_BPS_SECONDARY",method,.30)
                            complete=(all(results["modes"][c]["status"]=="NATIVE_MODES_VERIFIED"
                                         for c in ("slow","rapid")))
                            rec[method+"_three_layer"]=dict(
                                status="THREE_LAYER_NATIVE_TOV_F_P1_BOTH_CONVERSIONS_VERIFIED" if complete
                                       else "THREE_LAYER_NATIVE_TOV_PARTIAL_MODE_IDENTITY",
                                stellar=results,contract=str((dest/"REPRESENTATION.json").relative_to(ROOT)),
                                admission_note="A computed TOV is not full f/p1 validation if mode IDs are missing")
                        except Exception as ex:
                            rec[method+"_three_layer"]=dict(status="THREE_LAYER_BLOCKED",error=failure(ex))
                print("BBP_BPS_CASE",branch,name,rec["maxwell"]["status"],
                      rec["css"]["status"],
                      rec.get("maxwell_three_layer",{}).get("status"),
                      rec.get("css_three_layer",{}).get("status"),flush=True)
            except Exception as ex:
                rec=dict(branch=branch,core=name,status="SOURCE_FAILED",error=failure(ex))
            summary["combinations"].append(rec)
            save(OUT/"PARTIAL.json",summary)
    summary["status"]="EXECUTED_WITH_SOURCE_LEVEL_AND_STELLAR_GATES"
    save(OUT/"RESULT.json",summary)

if __name__=="__main__":main()
