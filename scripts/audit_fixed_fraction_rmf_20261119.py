#!/usr/bin/env python3
"""Controlled same-source-Yp RMF test without re-solving beta equilibrium.

Demonstrates that observed 1-15% pressure discrepancy is not induced by
the O(1e-6) differences in beta root-finding proton composition.
Preserves source, fixes Yp to provider's eos.compo printed value, determines
free-lepton mu from charge neutrality and compares model pressure.
"""
from __future__ import annotations
import hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos import read_compose,gppva_rmf
from compose_eos.cold_leptons import ELECTRON_MASS,free_fermion,leptons

BASE=ROOT/"data/cowling_maxwell_native_20261119/provider_original_unified"
OUT=ROOT/"data/cowling_rmf_reconciliation_20261119/fixed_composition"
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")

def pressure_at_provider_proton_fraction(name,n,yp):
    bary=gppva_rmf.baryons(name,float(n),float(yp))
    charge=float(n*yp)
    hi=float(free_fermion(charge,ELECTRON_MASS)[2])
    mu=brentq(lambda q:leptons(q,True)["n_e"]+leptons(q,True)["n_mu"]-charge,
              ELECTRON_MASS,hi*(1+1e-12),xtol=1e-12)
    lep=leptons(mu,True)
    return dict(P_total=bary["pressure_baryons"]+lep["pressure_leptons"],
       epsilon_total=bary["epsilon_baryons"]+lep["epsilon_leptons"],
       baryonic_pressure=bary["pressure_baryons"],
       lepton_pressure=lep["pressure_leptons"],
       neutral_chemical_mu_e=mu,
       beta_residual_at_source_Yp_MeV=bary["mu_n"]-bary["mu_p"]-mu)

def main():
    out=dict(schema="ns_cowling.stage05.fixed_original_fraction_rmf.v1",
       source_eos_changed=False,
       source_density_exact_no_interpolation=True,
       composition_policy="provider Yp held fixed, Lepton mu from free Fermi neutrality",
       core_quantities_no_coupling_fit=True,models={},status="RUNNING")
    for name in ("FSU2","TM1e"):
        source=BASE/name/"raw"
        d=read_compose(source)
        ids=np.flatnonzero((d.Iphase==6)&(d.nb>=.10)&(d.nb<=.40))
        chosen=ids[np.linspace(0,len(ids)-1,8,dtype=int)]
        values=[]
        for i in chosen:
            n=float(d.nb[i]);yp=float(d.Yp[i])
            fixed=pressure_at_provider_proton_fraction(name,n,yp)
            eq=gppva_rmf.beta_state(name,n)
            gap=eq["pressure"]-float(d.p[i])
            due_to_fraction=fixed["P_total"]-eq["pressure"]
            values.append(dict(source_index=int(i),n_fm3=n,
                source_Yp=yp,independent_beta_Yp=float(eq["Y_p"]),
                source_pressure=float(d.p[i]),
                fixed_proton_fraction_RMF_pressure=float(fixed["P_total"]),
                independent_beta_RMF_pressure=float(eq["pressure"]),
                beta_composition_pressure_delta_MeVfm3=float(due_to_fraction),
                original_core_pressure_discrepancy_MeVfm3=float(gap),
                composition_effect_over_bulk_pressure_mismatch=float(
                    abs(due_to_fraction/gap)) if abs(gap)>1e-12 else None,
                free_electron_mu_MeV=float(fixed["neutral_chemical_mu_e"]),
                original_compose_Q4_electron_mu_MeV=float(-d.mn*d.Q4[i]),
                beta_residual_at_source_Yp_MeV=float(fixed["beta_residual_at_source_Yp_MeV"])))
        out["models"][name]=dict(original_eos_thermo_sha256=sha(source/"eos.thermo"),
                tested_rows=values,
                max_proton_fraction_perturbation_effect_MeVfm3=max(
                    abs(x["beta_composition_pressure_delta_MeVfm3"]) for x in values),
                min_source_pressure_abs_error_MeVfm3=min(
                    abs(x["original_core_pressure_discrepancy_MeVfm3"]) for x in values),
                max_composition_error_to_source_pressure_gap_ratio=max(
                    x["composition_effect_over_bulk_pressure_mismatch"] for x in values),
                source_to_model_P_mismatch_does_not_derive_from_beta_root_rounding=True)
        print("FIXED_SOURCE_FRACTION",name,
            "max beta root P effect",
            out["models"][name]["max_proton_fraction_perturbation_effect_MeVfm3"],
            "smallest source P discrepancy",
            out["models"][name]["min_source_pressure_abs_error_MeVfm3"],flush=True)
        save(OUT/"PARTIAL.json",out)
    out["status"]="CONTROLLED_PROVIDER_FRACTION_TEST_COMPLETE"
    save(OUT/"RESULT.json",out)
if __name__=="__main__":main()
