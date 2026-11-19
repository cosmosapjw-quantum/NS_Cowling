#!/usr/bin/env python3
"""Source-locked original-GPPVA one-density and independent-author RMF audit.

E01: original Pais PDFs/ZIP provide no eos.micro, mean-field components or
     generator; source columns are compared exactly at original density.
E02: original Shen 2020 TM1e EOS4 author provides independent M*, P, F, mu
     on the exactly tabulated T=0, fixed Yp uniform grid (no leptons).
E03: local computed nonlinear field/P/epsilon decomposition, not falsely
     ascribed to original Pais generation, plus native DDME2 control.
E04: original FSU2 literature saturated symmetric M* and PNM pressure.

Do not modify EOS input, change parameters, or claim source-equivalence to
original Pais core generation where the original meson fields are unavailable.
"""
from __future__ import annotations
import csv,hashlib,json,sys,traceback
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos import read_compose,gppva_rmf
from compose_eos.native import evaluate_ddrmf
from compose_eos.cold_leptons import free_fermion,ELECTRON_MASS,MUON_MASS
from compose_eos.rmf_onepoint_components import (
    nonlinear_component_budget,read_shen_eos4_t00_exact_row)
SOURCE=ROOT/"data/cowling_maxwell_native_20261119/provider_original_unified"
OLD=ROOT/"data/cowling_rmf_reconciliation_20261119/EXACT_DENSITY_AUDIT.json"
SHEN=ROOT/"data/cowling_original_rmf_onepoint_20261119/external_reference/shen_EOS4_original_eos4.t00.gz"
OUT=ROOT/"data/cowling_original_rmf_onepoint_20261119"
TARGET_N=.2
FORMULAS=[
 "energy_total_baryons_and_free_leptons",
 "epsilon_vector = -W + omega*dW/domega + rho*dW/drho",
 "P_vector = W, scalar=-U, rearrangement=nB*Sigma_R",
 "sigma positive phi_local, sigma_author=-phi",
]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def annotate_source(name,source_hash):
    return dict(
        GPPVA_2023_CompOSE_Pais_source_url=f"https://compose.obspm.fr/eos/{'215' if name=='FSU2' else '221' if name=='TM1e' else '218'}",
        GPPVA_author="Helena Pais, University of Coimbra, hpais@uc.pt",
        original_GPPVA_eos_micro_available=False,
        original_GPPVA_software_generator_obtained=False,
        original_GPPVA_source_thermo_sha256=source_hash,
        source_status="ONE_EXACT_TABULATED_GPPVA_TOTAL_P_EPS_Y_AND_Q; NO_ORIGINAL_MESON_FIELDS")


def one(name,authority):
    orig=SOURCE/name/("validated_monotone_prefix_projection" if name=="DDME2" else "raw")
    base=SOURCE/name/"raw"
    authentic=sha(base/"eos.thermo")
    expected=authority[name]["source_identity"]["original_eos_thermo_sha256"]
    if authentic!=expected:
        raise ValueError("immutable original GPPVA source SHA drift: "+name)
    eos=read_compose(orig)
    inds=np.flatnonzero((eos.Iphase==6)&(eos.nb>=.08)&(eos.nb<=.45))
    row=int(inds[np.argmin(abs(eos.nb[inds]-TARGET_N))])
    n=float(eos.nb[row])
    if abs(n-TARGET_N)>.01:
        raise ValueError("no close original homogeneous source density")
    sp=dict(original_full_source_row_index_zero_based=row,
        original_source_n_B_fm3=n,
        original_pressure_MeV_fm3=float(eos.p[row]),
        original_epsilon_MeV_fm3=float(eos.eps[row]),
        original_mu_B_Q3_MeV=float(eos.mub[row]),
        original_mu_B_Euler_MeV=float((eos.p[row]+eos.eps[row])/n),
        original_charge_Q4_scaled=float(eos.Q4[row]),
        original_lepton_Q5_scaled=float(eos.Q5[row]),
        original_electron_mu_Q4Q5_MeV=float((eos.Q5[row]-eos.Q4[row])*eos.mn),
        original_Yp=float(eos.Yp[row]),original_Ye=float(eos.Ye[row]),
        original_Ymu=float(eos.Ymu[row]),
        original_particle_charge_residual_number_fraction=float(
             eos.Yp[row]-eos.Ye[row]-eos.Ymu[row]),
        original_neutron_reference_mass_MeV=float(eos.mn),
        source_muEuler_minus_muQ3_MeV=float(
            (eos.p[row]+eos.eps[row])/n-eos.mub[row]),
        original_electron_mu_from_Ye_MeV=float(
            free_fermion(n*eos.Ye[row],ELECTRON_MASS)[2]),
        original_electron_mu_Q_minus_Fermi_MeV=float(
             (eos.Q5[row]-eos.Q4[row])*eos.mn-
             free_fermion(n*eos.Ye[row],ELECTRON_MASS)[2]))
    if name=="DDME2":
        state=evaluate_ddrmf("GPPVA_DDME2",n)
        if state["status"]!="ok":
            raise RuntimeError("DDME2 native reference returned non-ok")
        local=dict(Yp=float(state["yp"]),Ye=float(state["ye"]),
            Ymu=float(state["ymu"]),
            pressure=float(state["pressure_total_mev_fm3"]),
            epsilon=float(state["eps_total_mev_fm3"]),
            muB=float(state["mu_n_mev"]),
            beta_residual=float(state["beta_residual_mev"]),
            identity_residual=float(state["pressure_identity_residual_mev_fm3"]),
            field_data_status="NOT_EXPOSED_BY_NATIVE_DDME2_REF_BACKEND",
            backend="existing native density-dependent Fortran RMF")
        components=None
    else:
        b=gppva_rmf.beta_state(name,n)
        local=dict(Yp=float(b["Y_p"]),Ye=float(b["Y_e"]),
            Ymu=float(b["Y_mu"]),pressure=float(b["pressure"]),
            epsilon=float(b["epsilon"]),muB=float(b["mu_B"]),
            beta_residual=float(b["beta_residual_mev"]),
            identity_residual=float(b["total_euler_residual"]),
            field_data_status="LOCAL_RMF_RECALCULATED_FIELDS_NOT_GPPVA_ORIGINAL",
            backend="independent stationary Python RMF")
        components=nonlinear_component_budget(name,n,b["Y_p"],
                                      n*b["Y_e"],n*b["Y_mu"])
        fixed=nonlinear_component_budget(name,n,float(eos.Yp[row]),
                         n*float(eos.Ye[row]),n*float(eos.Ymu[row]))
    compare=dict(
        provider_minus_local_pressure_MeV_fm3=float(sp["original_pressure_MeV_fm3"]-local["pressure"]),
        local_pressure_over_provider_minus1=float(local["pressure"]/sp["original_pressure_MeV_fm3"]-1),
        provider_minus_local_epsilon_MeV_fm3=float(sp["original_epsilon_MeV_fm3"]-local["epsilon"]),
        source_mub_Q3_minus_local_mu_B_MeV=float(sp["original_mu_B_Q3_MeV"]-local["muB"]),
        source_mub_Euler_minus_local_mu_B_MeV=float(sp["original_mu_B_Euler_MeV"]-local["muB"]),
        provider_minus_local_fraction=dict(Yp=float(sp["original_Yp"]-local["Yp"]),
             Ye=float(sp["original_Ye"]-local["Ye"]),
             Ymu=float(sp["original_Ymu"]-local["Ymu"])))
    result=dict(model=name,selected_density_is_original_provider_exact_knot=True,
         no_source_interpolation=True,source=sp,local_beta_equilibrium=local,
         source_minus_local=compare,source_authority=annotate_source(name,authentic),
         physical_gate="NORMALIZED_P_EPS_MU_Y_ONLY_NO_PAIS_MEAN_FIELDS")
    if components is not None:
        result["local_mean_field_components_at_local_beta_equilibrium"]=components
        result["local_mean_field_components_at_original_printed_Yp_Ye_Ymu"]=fixed
        result["source_minus_local_fixed_provider_Yp_pressure_MeV_fm3"]=float(
            sp["original_pressure_MeV_fm3"]-fixed["pressure_total_MeV_fm3"])
        result["source_minus_local_fixed_provider_Yp_epsilon_MeV_fm3"]=float(
            sp["original_epsilon_MeV_fm3"]-fixed["epsilon_total_MeV_fm3"])
        with (OUT/f"{name}_LOCAL_ONLY_COMPONENTS.csv").open("w",newline="") as f:
            writer=csv.writer(f)
            writer.writerow(["term","local_energy_MeV_fm3","local_pressure_MeV_fm3",
                             "original_GPPVA_component_status"])
            terms=sorted(set(components["epsilon_components_MeV_fm3"])|
                         set(components["pressure_components_MeV_fm3"]))
            for term in terms:
                writer.writerow([term,components["epsilon_components_MeV_fm3"].get(term,""),
                       components["pressure_components_MeV_fm3"].get(term,""),
                       "ORIGINAL_GPPVA_NOT_PUBLISHED"])
    return result


def original_shen_TM1e():
    r=read_shen_eos4_t00_exact_row(SHEN,14.5,.1)
    local=nonlinear_component_budget("TM1e",r["n_B_fm3"],r["Yp"],0.,0.)
    metrics=dict(
       Mstar_MeV=local["mstar_MeV"]-r["Mstar_MeV"],
       pressure_baryons_MeV_fm3=(local["pressure_baryons_MeV_fm3"]-
                                r["pressure_baryons_MeV_fm3"]),
       E_per_baryon_relative938_MeV=(local["epsilon_baryons_MeV_fm3"]/r["n_B_fm3"]-
                                 938.-r["free_energy_per_baryon_relative938_MeV"]),
       mu_n_relative938_MeV=local["mu_n_MeV"]-938.-r["mu_n_relative938_MeV"],
       mu_p_relative938_MeV=local["mu_p_MeV"]-938.-r["mu_p_relative938_MeV"],
       sigma_author_derived_from_Mstar_minus_local_MeV=(
          (r["Mstar_MeV"]-938.)/gppva_rmf.model("TM1e")["gs"]-
          local["local_sigma_author_sign_MeV"]))
    accepted=(abs(metrics["Mstar_MeV"])<.02
       and abs(metrics["pressure_baryons_MeV_fm3"])<.001
       and abs(metrics["E_per_baryon_relative938_MeV"])<.01
       and abs(metrics["mu_n_relative938_MeV"])<.05
       and abs(metrics["mu_p_relative938_MeV"])<.05)
    scalar_from_Mstar=(r["Mstar_MeV"]-938.)/gppva_rmf.model("TM1e")["gs"]
    return dict(original_author_source=r,local_fixed_Yp_baryons_only=local,
       author_implied_sigma_MeV=float(scalar_from_Mstar),
       local_sigma_in_author_sign_MeV=float(local["local_sigma_author_sign_MeV"]),
       derived_sigma_status="DERIVED_FROM_PUBLISHED_MSTAR_NOT_REPORTED_AS_FIELD_IN_ORIGINAL_TABLE",
       local_minus_shen_table=metrics,
       accepted_by_preregistered_tolerance=accepted,
       status=("INDEPENDENT_SHEN_TM1e_AUTHOR_BULK_EQUATION_VERIFIED"
               if accepted else "ORIGINAL_AUTHOR_RMF_COMPARISON_FAILED"),
       same_GPPVA_generator=False,
       caveat="Shen original functional not direct Pais GPPVA code or fields")


def FSU2_primary_reference():
    snm=gppva_rmf.baryons("FSU2",.1505,.5)
    pnm=gppva_rmf.baryons("FSU2",.1505,0.)
    return dict(
      source="Tolos, Centelles, Ramos, Astrophys J 834:3 (2017), Table 1",
      source_url="https://diposit.ub.edu/dspace/bitstream/2445/116038/1/669453.pdf",
      source_type="ORIGINAL_PUBLISHED_RMF_SATURATION_BENCHMARK_NOT_PAIS_GPPVA_ONEPOINT",
      n_B_fm3=.1505,
      original_FSU2_SNM_mstar_over_mN=.593,
      original_author_implied_scalar_phi_from_rounded_mstar_MeV=float(
           939*(1.-.593)/gppva_rmf.model("FSU2")["gs"]),
      local_scalar_phi_MeV=float(snm["phi"]),
      original_field_status="DERIVED_FROM_THREE_DIGIT_PUBLISHED_MSTAR_RATIO_NOT_AUTHOR_REPORTED_SIGMA",
      local_FSU2_SNM_mstar_over_mN=float(snm["effective_mass_n"]/939),
      original_FSU2_PNM_pressure_MeV_fm3=5.81,
      local_FSU2_PNM_pressure_MeV_fm3=float(pnm["pressure_baryons"]),
      status="SCOPED_LITERATURE_SATURATION_COMPARISON")

def main():
    authority=json.loads(OLD.read_text())["models"]
    OUT.mkdir(parents=True,exist_ok=True)
    report=dict(schema="ns_cowling.rmf_original_author_and_GPPVA_single_density.v1",
       experiment_ids=["E01","E02","E03","E04"],
       source_policy="original fixed density exact, no fit, no source repairs, no original field claim",
       source_05_checksum=sha(OLD),
       models={},independent_author_TM1e={},literature_FSU2={},
       normalized_fields=FORMULAS,original_fields_from_Pais_GPPVA="NOT_AVAILABLE",
       status="RUNNING")
    for name in ("FSU2","TM1e","DDME2"):
        report["models"][name]=one(name,authority)
        write(OUT/"PARTIAL.json",report)
        r=report["models"][name]
        print("GPPVA_EXACT_ONEPOINT",name,r["source"]["original_source_n_B_fm3"],
          "P_provider",r["source"]["original_pressure_MeV_fm3"],
          "P_local",r["local_beta_equilibrium"]["pressure"],
          "Pdiff",r["source_minus_local"]["provider_minus_local_pressure_MeV_fm3"],
          flush=True)
    report["independent_author_TM1e"]=original_shen_TM1e()
    report["literature_FSU2"]=FSU2_primary_reference()
    report["status"]="SCOPED_ONEPOINT_AUTHOR_COMPARISON_AND_LOCAL_TERMS_AUDITED"
    write(OUT/"RESULT.json",report)
    print("SHEN_ORIGINAL_AUTHOR_COMPARE",report["independent_author_TM1e"]["status"],
       report["independent_author_TM1e"]["local_minus_shen_table"],flush=True)
    print("SOURCE_ONEPOINT_RESULT",OUT/"RESULT.json",flush=True)

if __name__=="__main__":main()
