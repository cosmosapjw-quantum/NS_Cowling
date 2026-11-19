#!/usr/bin/env python3
"""Reproduce source calibration, including independent Zenodo byte evidence.

Compare the official homogeneous GPPVA phase6 rows at IDENTICAL n_B with
the local RMF function evaluator. No interpolation or EOS edits. Preserve the
difference between *same parameter family* and *source identical physics*.
"""
from __future__ import annotations
import hashlib,json,sys,traceback
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos import read_compose,gppva_rmf
from compose_eos.native import evaluate_ddrmf
from compose_eos.source_calibration import (first_law_segment_audit,
    source_model_row_audit,phase_contiguous_slices)

SOURCE=ROOT/"data/cowling_maxwell_native_20261119/provider_original_unified"
OLD=ROOT/"data/cowling_maxwell_native_20261119/provider_original_unified/ORIGINAL_CONTROLS.json"
OUT=ROOT/"data/cowling_eos_source_calibration_20261119/source_audit"
ZEN=ROOT/"data/cowling_eos_source_calibration_20261119/zenodo/ZENODO_INDEPENDENT_ARCHIVE_READBACK.json"
RECORDS=("FSU2","TM1e","DDME2")

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+"\n")

def actual_rmf_state(name,n):
    """Physically derived independent homogeneous EOS, no P(n) interpolation."""
    if name=="DDME2":
        x=evaluate_ddrmf("GPPVA_DDME2",n)
        if x["status"]!="ok":raise RuntimeError("native DDME2 RMF source state failed")
        return dict(n=n,P=float(x["pressure_total_mev_fm3"]),
            epsilon=float(x["eps_total_mev_fm3"]),mu_B=float(x["mu_n_mev"]),
            beta_residual_MeV=float(x["beta_residual_mev"]),
            firstlaw_pressure_identity_residual=float(x["pressure_identity_residual_mev_fm3"]),
            evaluator="Fortran DD-RMF")
    x=gppva_rmf.beta_state(name,n)
    if abs(x["total_euler_residual"])>1e-7:
        raise ValueError("Python RMF homogeneous reference violates Euler")
    return dict(n=n,P=float(x["pressure"]),epsilon=float(x["epsilon"]),
        mu_B=float(x["mu_B"]), beta_residual_MeV=float(x["beta_residual_mev"]),
        firstlaw_pressure_identity_residual=float(x["total_euler_residual"]),
        evaluator="Python original published-RMF parameter solver")

def compact_reconstruction(a):
    return {k:v for k,v in a.items()
            if k not in ("n_source","p_source","epsilon_source","source_mu_Q3",
                         "source_mu_Euler","corrected_pressure","corrected_epsilon",
                         "corrected_mu","corrected_f_energy_per_baryon")}

def audit(name,original,zen,archive):
    provider=SOURCE/name/"raw"
    projection=(ROOT/original["source_projection"]["projection_path"]
                if original.get("source_projection") else provider)
    d=read_compose(projection)
    fixed_zen=zen[name]
    observed_zip=provider.parent/"eos.zip"
    if digest(observed_zip)!=fixed_zen["observed_zip_sha256"]:
        raise RuntimeError("original CompOSE bytes changed vs Zenodo download")
    if digest(provider/"eos.thermo")!=original["source_file_hashes"]["eos.thermo"]:
        raise RuntimeError("official source thermo bytes changed")
    phasegroups=phase_contiguous_slices(d.Iphase)
    segments=[]
    for phase,start,stop in phasegroups:
        group=dict(Iphase=int(phase),first_row=int(start),stop_row_exclusive=int(stop),
            density_support=[float(d.nb[start]),float(d.nb[stop-1])],
            pressure_support=[float(d.p[start]),float(d.p[stop-1])])
        if stop-start>=3:
            try:
                audit=first_law_segment_audit(d.nb[start:stop],d.p[start:stop],
                                              d.eps[start:stop],d.mub[start:stop])
                group["first_law"]=compact_reconstruction(audit)
            except Exception as ex:
                group["first_law"]=dict(status="BLOCKED_SOURCE_LOCAL",
                      type=type(ex).__name__,error=str(ex)[:300])
        segments.append(group)
    coreidx=np.flatnonzero((np.asarray(d.Iphase)==6)
                            &(d.nb>=0.08)&(d.nb<=.45))
    if len(coreidx)<15:raise RuntimeError("too few official core points in calibration support")
    evaluator=[actual_rmf_state(name,float(d.nb[i])) for i in coreidx]
    src=dict(n=d.nb[coreidx],P=d.p[coreidx],
             epsilon=d.eps[coreidx],mu_Q3=d.mub[coreidx])
    cmp=source_model_row_audit(src,evaluator)
    # Percent-level differences must not be labelled source-identical merely
    # because the same nuclear functional name appears in both pathways.
    cmp["classification"]=("HIGH_NUMERIC_AGREEMENT_NOT_BIT_IDENTITY"
          if cmp["pressure_abs_rel_max"]<.0002 and cmp["epsilon_abs_rel_max"]<.0002
          else "SAME_NAMED_PARAMETER_FAMILY_BUT_SOURCE_NOT_NUMERICALLY_EQUIVALENT")
    cmp["evaluator"]=evaluator[0]["evaluator"]
    cmp["max_abs_beta_residual_MeV"]=max(abs(x["beta_residual_MeV"]) for x in evaluator)
    cmp["max_abs_euler_residual_MeV_fm3"]=max(abs(x["firstlaw_pressure_identity_residual"]) for x in evaluator)
    return dict(model=name,original_homogeneous_phase_count=sum(x[0]==6 for x in phasegroups),
        provider_raw_source_full_n=int(original.get("original_layout",{}).get("total_rows",
                                   original.get("n_rows_source",d.n))),
        analyzed_n=int(d.n),
        source_trim_explanation=original.get("source_projection"),
        zipped_source_archive_sha256=digest(observed_zip),
        independent_zenodo_byte_identity=True,
        independent_zenodo_record_url="https://zenodo.org/records/"+str(fixed_zen["record"]),
        provider_listing_sha256_disagreement_retained=True,
        raw_eos_source_sha256={f.name:digest(f) for f in provider.glob("eos.*") if f.is_file()},
        phase_contiguous_audits=segments,
        original_homogeneous_core_same_density_comparison=cmp)

def main():
    zen=json.loads(ZEN.read_text());old=json.loads(OLD.read_text())
    out=dict(schema="ns_cowling.eos_source_calibration.original_vs_same_named_rmf.v1",
       units=dict(n="fm^-3",P="MeV/fm^3",epsilon="MeV/fm^3",mu_B="MeV"),
       no_changes_to_original_EOS=True,
       original_Q3_Euler_separated_from_first_law_reconstruction=True,
       provider_archive_secondary_authentication="Zenodo independent zipped bytes",
       archive_receipt_path=str(ZEN.relative_to(ROOT)),
       sample_density_range_fm3=[0.08,0.45],
       original_source_model_equivalence_not_assumed=True,
       models={},status="RUNNING")
    for name in RECORDS:
        try:
            data=audit(name,old["models"][name],zen,old)
            out["models"][name]=data
            print("SOURCE_CALIBRATION",name,
                "pressure_max_relative",data["original_homogeneous_core_same_density_comparison"]["pressure_abs_rel_max"],
                "energy_max_relative",data["original_homogeneous_core_same_density_comparison"]["epsilon_abs_rel_max"],
                "mu_Q3_gap_max",data["original_homogeneous_core_same_density_comparison"]["mu_model_minus_provider_q3_abs_max_MeV"],
                "phase_count",len(data["phase_contiguous_audits"]),flush=True)
        except Exception as ex:
            out["models"][name]=dict(status="BLOCKED",error=str(ex),
                type=type(ex).__name__,traceback=traceback.format_exc()[-1200:])
            print("SOURCE_CALIBRATION_BLOCKED",name,type(ex).__name__,str(ex),flush=True)
        save(OUT/"PARTIAL.json",out)
    out["status"]="ALL_AVAILABLE_SOURCES_AUDITED"
    save(OUT/"RESULT.json",out)
    print("SOURCE_CALIBRATION_DONE",len(out["models"]),flush=True)

if __name__=="__main__":main()
