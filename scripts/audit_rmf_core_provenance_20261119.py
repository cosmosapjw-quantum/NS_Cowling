#!/usr/bin/env python3
"""Stage-05 exact-density core reconciliation and discriminating controls.

Original CompOSE FSU2, TM1e, DDME2 EOS source bytes immutable.
Examine beta-equilibrium fractions, p(epsilon) thermodynamics, symmetric
saturation, original Q3 vs Euler, and one literature-supplied precision
variant for TM1e, without parameter fitting or causal claim inflation.
"""
from __future__ import annotations
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos import read_compose,gppva_rmf
from compose_eos.native import evaluate_ddrmf
from compose_eos.rmf_source_reconciliation import (
   compare_phase_homogeneous_states,cold_pressure_from_energy,
   rest_mass_shift_pressure_invariant,decompose_pressure_difference,
)

BASE=ROOT/"data/cowling_maxwell_native_20261119/provider_original_unified"
OUT=ROOT/"data/cowling_rmf_reconciliation_20261119"
STAGE04=ROOT/"data/cowling_eos_source_calibration_20261119/source_audit/RESULT.json"
ZENODO=ROOT/"data/cowling_eos_source_calibration_20261119/zenodo/ZENODO_INDEPENDENT_ARCHIVE_READBACK.json"
TM1_PRECISE=dict(ms=511.19777,gs=10.02892,gw=12.61394,g2=-7.23247,g3=.61833,c3=71.30747)
REFERENCE={
    "FSU2":dict(n0=.1505,E_bind_MeV=-16.28,K_MeV=238.,source_url="https://compose.obspm.fr/eos/215"),
    "TM1e":dict(n0=.145,E_bind_MeV=-16.3,K_MeV=281.,source_url="https://compose.obspm.fr/eos/221"),
    "DDME2":dict(n0=.152,E_bind_MeV=-16.14,K_MeV=251.,source_url="https://compose.obspm.fr/eos/218"),
}

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n")

def point(name,n):
    if name=="DDME2":
        x=evaluate_ddrmf("GPPVA_DDME2",float(n))
        if x["status"]!="ok":raise RuntimeError("native RMF DDME2 evaluated state rejected")
        return dict(n=float(n),P=float(x["pressure_total_mev_fm3"]),
            epsilon=float(x["eps_total_mev_fm3"]),mu_B=float(x["mu_n_mev"]),
            Yp=float(x["yp"]),Ye=float(x["ye"]),Ymu=float(x["ymu"]),
            beta_residual_MeV=float(x["beta_residual_mev"]),
            euler_residual_MeV_fm3=float(x["pressure_identity_residual_mev_fm3"]))
    x=gppva_rmf.beta_state(name,float(n))
    return dict(n=float(n),P=float(x["pressure"]),epsilon=float(x["epsilon"]),
            mu_B=float(x["mu_B"]),Yp=float(x["Y_p"]),Ye=float(x["Y_e"]),
            Ymu=float(x["Y_mu"]), beta_residual_MeV=float(x["beta_residual_mev"]),
            euler_residual_MeV_fm3=float(x["total_euler_residual"]))

def saturation(name):
    n0=brentq(lambda n:gppva_rmf.baryons(name,float(n),.5)["pressure_baryons"],.08,.22,
              xtol=2e-13)
    state=gppva_rmf.baryons(name,n0,.5)
    step=2e-4
    grad=(-gppva_rmf.baryons(name,n0+2*step,.5)["pressure_baryons"]
        +8*gppva_rmf.baryons(name,n0+step,.5)["pressure_baryons"]
        -8*gppva_rmf.baryons(name,n0-step,.5)["pressure_baryons"]
        +gppva_rmf.baryons(name,n0-2*step,.5)["pressure_baryons"])/(12*step)
    return dict(n0_fm3=float(n0),E_bind_MeV=float(state["epsilon_baryons"]/n0-
                                  gppva_rmf.model(name)["mn"]),
                K_MeV=float(9*grad),effective_mass_over_rest=float(
                      state["effective_mass_n"]/gppva_rmf.model(name)["mn"]),
                published_saturation=REFERENCE[name],
                compared_convention="symmetric cold matter Yp=0.5 without leptons")

def source_original(name):
    base=BASE/name
    original=base/("validated_monotone_prefix_projection" if name=="DDME2" else "raw")
    source=read_compose(original)
    authoritative_hash=json.loads(STAGE04.read_text())["models"][name][
        "raw_eos_source_sha256"]["eos.thermo"]
    if digest(base/"raw/eos.thermo")!=authoritative_hash:
        raise RuntimeError("original source hash mismatch to published prior calibration")
    return source,dict(path=str(original.relative_to(ROOT)),
               original_eos_thermo_sha256=authoritative_hash,
               conditional_truncation_used=(name=="DDME2"))

def audit(name):
    t=time.monotonic()
    src,identity=source_original(name)
    phase=np.asarray(src.Iphase)
    indexes=np.flatnonzero((phase==6)&(src.nb>=.08)&(src.nb<=.45))
    if len(indexes)<40:raise RuntimeError("insufficient original core coverage")
    n=src.nb[indexes]
    models=[point(name,float(nv)) for nv in n]
    source=dict(n=n,P=src.p[indexes],epsilon=src.eps[indexes],
               mu_Q3=src.mub[indexes],Yp=src.Yp[indexes],
               Ye=src.Ye[indexes],Ymu=src.Ymu[indexes])
    comparison=compare_phase_homogeneous_states(source,models)
    thermo=cold_pressure_from_energy(n,source["epsilon"],source["P"],source["mu_Q3"])
    local_derivative=cold_pressure_from_energy(n,
         [m["epsilon"] for m in models],[m["P"] for m in models],
         [m["mu_B"] for m in models])
    ctrl=rest_mass_shift_pressure_invariant(n,source["epsilon"],source["P"],
                                           delta_mass_MeV=1.0)
    muon_control=None
    if name!="DDME2":
        idx=-1
        muon_free=gppva_rmf.beta_state(name,float(n[idx]),include_muons=False)
        muon_control=dict(n_fm3=float(n[idx]),
            provider_Ymu=float(source["Ymu"][idx]),
            with_muons_Ymu=float(models[idx]["Ymu"]),
            forced_no_muon_Yp=float(muon_free["Y_p"]),
            provider_Yp=float(source["Yp"][idx]),
            forced_no_muon_Yp_gap=float(muon_free["Y_p"]-source["Yp"][idx]),
            no_muon_pressure_diff_MeV_fm3=float(
                muon_free["pressure"]-source["P"][idx]),
            source_composition_has_muons=bool(source["Ymu"][idx]>0))
    precision=None
    if name=="TM1e":
        keep=gppva_rmf._NL["TM1e"].copy()
        try:
            gppva_rmf._NL["TM1e"].update(TM1_PRECISE)
            selected=np.asarray([0,len(n)//2,len(n)-1])
            responses=[point("TM1e",float(n[i])) for i in selected]
            precision=dict(source_precision_reference_url=
                   "https://user.numazu-ct.ac.jp/~sumi/eos/table2/guide_EOS1.pdf",
              high_precision_isoscalar_parameter_fields=TM1_PRECISE,
              no_parameter_fitting=True,
              densities=[float(n[i]) for i in selected],
              baseline_pressure_error=[float(models[i]["P"]-source["P"][i]) for i in selected],
              high_precision_pressure_error=[float(row["P"]-source["P"][i]) for i,row in zip(selected,responses)],
              baseline_energy_error=[float(models[i]["epsilon"]-source["epsilon"][i]) for i in selected],
              high_precision_energy_error=[float(row["epsilon"]-source["epsilon"][i]) for i,row in zip(selected,responses)])
        finally:
            gppva_rmf._NL["TM1e"].clear()
            gppva_rmf._NL["TM1e"].update(keep)
    split=decompose_pressure_difference(n,source["P"],source["epsilon"],
                 [m["P"] for m in models],[m["epsilon"] for m in models])
    split_summary={k:v for k,v in split.items()
                   if isinstance(v,(str,bool,int,float))}
    split_summary["sampled_source_density_rows"]=[
       dict(n_fm3=float(n[i]),
        observed_P_difference=float(split["pressure_observed_difference"][i]),
        source_energy_functional_derivative=float(split["pressure_from_energy_difference"][i]),
        source_firstlaw_identity_residual=float(split["provider_identity_residual"][i]),
        local_RMF_firstlaw_identity_residual=float(split["model_identity_residual"][i]))
       for i in range(len(n))]
    median_pressure=float(np.median(np.asarray([v["pressure_delta"] for v in comparison["homogeneous_rows"]])))
    out=dict(model=name,source_identity=identity,
        original_compose_archive_authentic_record=json.loads(ZENODO.read_text())[name],
        parameters=gppva_rmf.model(name),
        original_beta_equilibrium_fraction_comparison=comparison,
        provider_thermodynamic_derivative=dict(
            max_inner_Pminus_n2dF_abs_MeV_fm3=thermo["max_inner_abs_P_error"],
            median_signed_Pminus_n2dF_MeV_fm3=thermo["median_signed_inner_PminusPderived_MeV_fm3"],
            median_abs_rel_Pminus_n2dF=thermo["median_inner_abs_rel_P_error"],
            max_abs_Q3minus_Euler_MeV=thermo["max_provider_Q3_Euler_discrepancy_MeV"],
            caution="cubic sampled derivative, exclude two endpoints: not a certified repaired EOS"),
        independent_local_thermodynamic_derivative=dict(
            median_abs_rel_Pminus_n2dF=local_derivative["median_inner_abs_rel_P_error"],
            max_abs_Q3minus_Euler_MeV=local_derivative["max_provider_Q3_Euler_discrepancy_MeV"]),
        pressure_difference_decomposition=split_summary,
        only_constant_baryon_mass_offset_invalid=dict(
            pressure_invariant=ctrl["max_abs_pressure_shift"]==0,
            actual_bulk_P_mismatch_max_relative=comparison["pressure_max_abs_rel"],
            explanation="epsilon->epsilon+n*delta_m changes Euler mu by constant but leaves pressure exactly invariant"),
        positive_median_model_minus_provider_pressure=bool(median_pressure>0),
        symmetric_matter_saturation=saturation(name),
        lepton_species_hypothesis=muon_control,
        parameter_table_precision_hypothesis=precision,
        experiment_seconds=time.monotonic()-t,
        source_equivalence_status="NOT_CONFIRMED_NONZERO_PRESSURE_DIFFERENCE",
        no_fitted_coupling_or_energy_offset=True)
    print("AUDIT",name,"n",len(n),
      "max Y",comparison["max_abs_composition_fraction_difference"],
      "rounding match",comparison["composition_compatible_with_printed_rounding"],
      "max P%",100*comparison["pressure_max_abs_rel"],
      "relative P-from-eps median",out["provider_thermodynamic_derivative"][
            "median_abs_rel_Pminus_n2dF"],flush=True)
    return out

def main():
    output=dict(schema="ns_cowling.rmf_source_reconciliation_05.full_phase_at_identical_density.v1",
       material_claim="compositional source match does not prove full source-identical bulk EOS",
       method="independent exact source n RMF beta solutions, source EOS immutable",
       source_04_identity_sha256=digest(STAGE04),
       source_04_Zenodo_hash_receipt_sha256=digest(ZENODO),
       models={},status="RUNNING")
    for model_name in ("FSU2","TM1e","DDME2"):
        output["models"][model_name]=audit(model_name)
        save(OUT/"PARTIAL.json",output)
    output["status"]="BOUNDED_ORIGINAL_SOURCE_RECONCILIATION_AUDITED"
    save(OUT/"EXACT_DENSITY_AUDIT.json",output)
    print("DONE",[(k,v["source_equivalence_status"]) for k,v in output["models"].items()],flush=True)

if __name__=="__main__":main()
