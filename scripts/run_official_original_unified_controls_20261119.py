#!/usr/bin/env python3
"""Official GPPVA unified EOS readback, phase split, native Maxwell/CSS
attempts and independent minimally thermodynamically corrected smooth star.

Do not promote a failure of an original table to a physics no-go.
No uncontrolled extrapolation across the chemical potential support gap.
"""
import hashlib
import json
import sys
import traceback
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos import read_compose,match_crust_core
from compose_eos.api import construct_maxwell_transition
from compose_eos.maxwell_v2 import _segment_eosdata
from compose_eos.whole_crust import (
    SurfacePolytrope,PiecewisePolytrope,CausalExponentialInterval,WholeCrustBarotrope)
from compose_eos.stellar_background import ContractEnthalpyEOS,solve_background
from compose_eos.stellar_cowling import CowlingSettings
from compose_eos.stellar_cowling_native import NativeCowlingDVODE

ROOT_DATA=ROOT/"data/cowling_maxwell_native_20261119"
OUT=ROOT_DATA/"provider_original_unified"
MODELS=("FSU2","TM1e","DDME2")
N_C=.30
N_MAX=.49

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,value):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def error(exc):return dict(type=type(exc).__name__,message=str(exc)[:520])

def _view(source, start, stop, *, first_law=None):
    """Explicit purely thermodynamic original interval projection.

    Preserve n and P pointwise and original ε and μ unless a documented,
    separately saved pressure-primary first-law reconstruction is requested.
    """
    obj=_segment_eosdata(source,(0,int(start),int(stop)))
    obj.composition_available=False;obj.composition_file_present=False
    obj.q_extra_count=0
    obj.q_extra.fill(0)
    obj.Iphase[:]=1
    for name in ("Yn","Yp","Ye","Ymu","Aav","Zav","Yh"):
        getattr(obj,name)[:]=0
    for name in ("Yn_present","Yp_present","Ye_present","Ymu_present",
                 "Abar_present","Zbar_present","Yh_present","C3_applicable",
                 "C4_applicable"):
        getattr(obj,name)[:]=False
    obj.row_status[:]=0
    if first_law is not None:
        obj.Q3[:]=first_law.mu[obj.nb.searchsorted(obj.nb)].copy()/obj.mn-1.
        obj.Q6[:]=first_law.F[obj.nb.searchsorted(obj.nb)].copy()/obj.mn-1.
        obj.Q7[:]=obj.Q6.copy()
    return obj

def attempt(l,h,method):
    try:
        if method=="maxwell":
            out,rep=construct_maxwell_transition(l,h,
                 assume_local_charge_neutrality=True,drop_composition=True,
                 return_report=True)
            return dict(status="PASS_NATIVE_FINITE",
                        transition_pressure=rep["pressure_transition"],
                        transition_mu=rep["mu_transition"],
                        density_low=rep["n_low"],density_high=rep["n_high"])
        out,rep=match_crust_core(l,h,"ORIGINAL_INNER","ORIGINAL_CORE",
                                "css",40,drop_composition=True)
        return dict(status="PASS_NATIVE_FINITE",
                    transition_pressure=rep["css_transition_pressure"],
                    density_low=rep["css_transition_nb_low"],
                    density_high=rep["css_transition_nb_high"])
    except Exception as exc:
        return dict(status="NATIVE_REJECTED",failure=error(exc))

def run_model(model):
    folder=OUT/model
    archive=folder/"eos.zip"
    checksum=folder/"eos.zip_checksum.txt"
    raw=folder/"raw"
    row=dict(model=model,
        official_source_url=f"https://compose.obspm.fr/download/1D/NS/RMF/{model}/eos.zip",
        provider_eos_zip_sha256=sha(archive),
        provider_checksum_recorded=checksum.read_text().split()[0].strip(),
        provider_zip_checksum_matches=(sha(archive)==checksum.read_text().split()[0].strip()),
        source_file_hashes={p.name:sha(p) for p in raw.glob("eos.*") if p.is_file()},
        comparison_central_density_fm3=N_C, comparison_upper_density_cutoff_fm3=N_MAX,
        states={})
    independent=folder/"eos.thermo.independent_download"
    row["independent_thermo_download_sha256"]=sha(independent) if independent.exists() else None
    row["independent_thermo_matches_zip_member"]=(sha(independent)==sha(raw/"eos.thermo")) if independent.exists() else None
    try:
        data=read_compose(raw)
    except Exception as exc:
        nb=(raw/"eos.nb").read_text().splitlines()
        densities=np.asarray([float(s) for s in nb[2:] if s.strip()])
        violations=np.flatnonzero(np.diff(densities)<=0)
        row.update(source_full_read_failure=error(exc),
                   n_rows_source=int(len(densities)),
                   density_nonmonotone_source_indices=violations.tolist(),
                   exact_original_not_overwritten=True)
        # Only a strictly monotone original prefix can be projected. Never
        # reorder, spline, extrapolate, or silently discard interior rows.
        if model!='DDME2' or len(violations)==0 or violations[0]<10:
            row['status']='RAW_FULL_SOURCE_BLOCKED';return row
        stop=int(violations[0]+1)
        dest=folder/'validated_monotone_prefix_projection'
        dest.mkdir(parents=True,exist_ok=True)
        for suffix in ('nb','thermo','compo'):
            lines=(raw/('eos.'+suffix)).read_text().splitlines()
            if suffix=='nb':
                assert int(lines[1].strip())==len(densities)
                replacement=[lines[0],str(stop),*lines[2:2+stop]]
            elif suffix=='thermo':
                assert len(lines)==len(densities)+1
                replacement=lines[:stop+1]
            else:
                assert len(lines)==len(densities)
                replacement=lines[:stop]
            (dest/('eos.'+suffix)).write_text('\n'.join(replacement)+'\n')
        for suffix in ('t','yq'):
            import shutil
            shutil.copyfile(raw/('eos.'+suffix),dest/('eos.'+suffix))
        data=read_compose(dest)
        if data.n!=stop or not np.array_equal(data.nb,densities[:stop]):
            raise RuntimeError('monotone prefix projection changes original density bytes')
        row['source_projection']=dict(
            kind='VALIDATED_UNCHANGED_PREFIX_NOT_ORIGINAL_COMPLETE',
            retained_source_rows=stop,original_rows=len(densities),
            rejected_tail_source_rows=list(range(stop+1,len(densities)+1)),
            physical_density_support_up_to_fm3=float(densities[stop-1]),
            projection_path=str(dest.relative_to(ROOT)),
            projected_eos_thermo_sha256=sha(dest/'eos.thermo'),
            truncation_reason='density inversion in last source rows: no extrapolation or sorting',
            raw_eos_thermo_sha256=sha(raw/'eos.thermo'))
        print('OFFICIAL_DDME2_SOURCE_PREFIX',stop,len(densities),flush=True)
    ph=np.asarray(data.Iphase)
    inds0=np.flatnonzero(ph==0)
    inds7=np.flatnonzero(ph==7)
    inds6=np.flatnonzero(ph==6)
    if not (len(inds0)>5 and len(inds7)>5 and len(inds6)>5 and
            inds0[-1]+1==inds7[0] and inds7[-1]+1==inds6[0]):
        row["status"]="SOURCE_LAYER_LAYOUT_UNRESOLVED"
        return row
    i=int(inds0[-1]+1);j=int(inds6[0])
    row["original_layout"]=dict(total_rows=int(data.n),outer_rows=i,
        inner_rows=j-i,core_rows=int(data.n-j),
        original_seams=[dict(left_n=float(data.nb[k-1]),right_n=float(data.nb[k]),
                              left_p=float(data.p[k-1]),right_p=float(data.p[k]),
                              left_mu=float(data.mub[k-1]),right_mu=float(data.mub[k])) for k in (i,j)])
    lo_raw=_view(data,i,j)
    hi_raw=_view(data,j,data.n)
    row["original_split_mu_overlap"]=[
        float(max(lo_raw.mub[0],hi_raw.mub[0])),
        float(min(lo_raw.mub[-1],hi_raw.mub[-1]))]
    row["original_split_has_mu_overlap"]=bool(row["original_split_mu_overlap"][1]>
                                               row["original_split_mu_overlap"][0])
    row["raw_original_split_methods"]={
        method:attempt(lo_raw,hi_raw,method) for method in ("maxwell","css")}
    cutoff=int(np.searchsorted(data.nb,N_MAX,side="right"))
    if cutoff<=j+6:
        row["status"]="NO_SUFFICIENT_TRUNCATED_CORE_SUPPORT"
        return row
    # Minimal thermodynamic correction: only F=ε/n and μ=F+P/n follow from
    # dF/dn=P/n². EVERY SOURCE n,P IS UNCHANGED. This changes physical ε
    # and μ, so it is never represented as the untouched original EOS.
    outer=PiecewisePolytrope(data.nb[:i],data.p[:i],data.eps[0]/data.nb[0])
    upper=PiecewisePolytrope(data.nb[i:cutoff],data.p[i:cutoff],data.eps[i]/data.nb[i])
    if outer.validate()["status"]!="PASS" or upper.validate()["status"]!="PASS":
        row["status"]="PRESSURE_PRIMARY_FIRST_LAW_FAILED"
        return row
    row["first_law_corrections"]=dict(
       outer_max_abs_dF_MeV=float(np.max(abs(outer.F-data.eps[:i]/data.nb[:i]))),
       upper_max_abs_dF_MeV=float(np.max(abs(upper.F-data.eps[i:cutoff]/data.nb[i:cutoff]))),
       upper_max_abs_relative_de=float(np.max(abs((upper.F*data.nb[i:cutoff]-data.eps[i:cutoff])
                                                /data.eps[i:cutoff]))),
       density_pressure_source_knots_preserved=True,
       anchor_each_region_original_first_eps_over_n=True,
       source_upper_truncated_to_support_n_leq=N_MAX,
       global_energy_offset_fitted=False)
    corrected=data.copy()
    corrected.Q3[i:cutoff]=upper.mu/corrected.mn-1
    corrected.Q6[i:cutoff]=upper.F/corrected.mn-1
    corrected.Q7[i:cutoff]=corrected.Q6[i:cutoff]
    l=_view(corrected,i,j)
    h=_view(corrected,j,cutoff)
    row["first_law_corrected_split_methods"]={
        method:attempt(l,h,method) for method in ("maxwell","css")}
    surface=SurfacePolytrope.from_anchor(data.nb[0],data.p[0],data.eps[0])
    bridge=CausalExponentialInterval(float(outer.nb[-1]),float(outer.p[-1]),
          float(outer.nb[-1]*outer.F[-1]),float(data.nb[i]),float(data.p[i]),
          float(data.eps[i]))
    bvalidation=bridge.validate()
    if bvalidation["status"]!="PASS":
        row["status"]="LOW_DENSITY_CONNECTOR_REJECTED"
        return row
    original_raw_upper=dict(n=data.nb[i:cutoff].tolist(),
         p=data.p[i:cutoff].tolist(),epsilon=data.eps[i:cutoff].tolist(),
         sharp_interfaces=[])
    record=dict(surface=surface.contract(),outer=outer.contract(),
           lower_bridge=bridge.contract(),upper=original_raw_upper)
    try:
        _=WholeCrustBarotrope.from_contract(record)
        row["raw_hermite_whole_star_status"]="PASS"
    except Exception as exc:
        row["raw_hermite_whole_star_status"]="REJECTED:"+str(exc)[:130]
    record["upper"]["epsilon"]=(upper.F*data.nb[i:cutoff]).tolist()
    whole=WholeCrustBarotrope.from_contract(record)
    target=folder/"thermodynamically_corrected_smooth"
    save(target/"REPRESENTATION.json",record)
    eos=ContractEnthalpyEOS(whole)
    control=dict(surface_h=1e-7,rtol=2e-9,atol=2e-11,max_step_h=.004)
    bg=solve_background(eos,N_C,backend="fortran_dvode",native_samples=8193,
                         **control)
    modes=NativeCowlingDVODE(bg,CowlingSettings(surface_h=1e-7,rtol=3e-9,
                atol=3e-12,max_step_h=.001),samples=8193)
    brackets,_=modes.scan_brackets(np.linspace(1000,9000,81))
    roots=[]
    for br in brackets[:6]:
        out=modes.solve_mode(br,profile_points=260)
        if out["W_nodes"] in (0,1) and out["V_nodes"]==out["W_nodes"]:
            roots.append(dict(label="f" if out["W_nodes"]==0 else "p1",
                frequency_hz=out["frequency_hz"],W_nodes=out["W_nodes"],
                V_nodes=out["V_nodes"],matching_residual=out["matching_residual"],
                bracket_hz=list(br)))
        if len({x["label"] for x in roots})>=2:break
    row["corrected_smooth_native"]=dict(status="PASS" if len(roots)>=2 else "PARTIAL",
         source_is_original_unified_pressure_primary_first_law_corrected=True,
         not_exact_raw_thermodynamic_eos=True,
         lower_gap_has_explicit_causal_connector=True,
         center_density_fm3=N_C,mass_msun=bg.mass_msun,radius_km=bg.radius_km,
         modes=roots,report_contract=str((target/"REPRESENTATION.json").relative_to(ROOT)),
         first_law_bridge_validation=bvalidation)
    if model=="FSU2":
        matched=json.loads((ROOT_DATA/"three_layer_original_fsu2/RESULT.json").read_text())
        comparisons={}
        bylabel={r["label"]:r["frequency_hz"] for r in roots}
        for method,decision in matched["methods"].items():
            star=decision.get("admitted",{}).get("stellar")
            if not star:continue
            cases={}
            for conv in ("slow","rapid"):
                mode_dict={z["label"]:z["frequency_hz"] for z in star["modes"][conv]["modes"]}
                cases[conv]={name:mode_dict.get(name,float("nan"))-bylabel.get(name,float("nan"))
                             for name in ("f","p1")}
            comparisons[method]=dict(delta_M_msun=star["native_tov"]["mass_msun"]-bg.mass_msun,
                delta_R_km=star["native_tov"]["R_km"]-bg.radius_km,
                delta_frequency_hz=cases,
                effects_not_isolated="different original source branch, thermodynamic corrections, finite Maxwell/CSS jump and connector")
        row["difference_from_same_rmf_reconstructed_sharp"]=comparisons
    row["status"]="ORIGINAL_UNIFIED_SPLIT_AND_CORRECTED_STELLAR_CONTROL_VERIFIED"
    return row

def main():
    output=dict(schema="ns_cowling.official_unified_original_compare.v1",
      status="RUNNING",scope="bounded source-original split + optional first-law corrected controls",
      nuclear_matter="original GPPVA full-source EOS from official provider zip",
      source_checksum_policy="Official ZIP checksum mismatch is never silently ignored; member EOS bytes independently checked when available",
      models={})
    for name in MODELS:
        try:
            r=run_model(name)
            output["models"][name]=r
            print("OFFICIAL_ORIGINAL",name,r["status"],
                  r.get("raw_original_split_methods"),
                  r.get("first_law_corrected_split_methods"),
                  r.get("corrected_smooth_native",{}).get("modes"),flush=True)
        except Exception as exc:
            output["models"][name]=dict(status="UNRESOLVED_IMPLEMENTATION",
                failure=error(exc),traceback=traceback.format_exc()[-1800:])
            print("ORIGINAL_UNRESOLVED",name,type(exc).__name__,str(exc),flush=True)
        save(OUT/"ORIGINAL_CONTROLS_PARTIAL.json",output)
    output["status"]="BOUNDED_EXPERIMENT_RECORDED"
    save(OUT/"ORIGINAL_CONTROLS.json",output)
    print("OFFICIAL_ORIGINAL_CONTROLS_FINAL",len(output["models"]),flush=True)

if __name__=="__main__":main()
