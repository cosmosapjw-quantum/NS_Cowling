#!/usr/bin/env python3
"""Generate homogeneous GDFM candidate CompOSE tables with explicit npe-muon composition."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.api import check_consistency, read_compose, write_compose
from compose_eos.eos_data import EOSData
from compose_eos.gdfm_uniform import MASS, MODEL_IDS, baryons, beta_state, parameters, symmetry_energy
from compose_eos.physical_domain import physical_domain_report


def write_model(model: str, density: np.ndarray, output: Path) -> dict:
    output.mkdir(parents=True,exist_ok=False)
    receipt = {"model_id":model,"status":"RUNNING"}
    receipt_path = output/"generation_receipt.json"
    receipt_path.write_text(json.dumps(receipt)+"\n")
    try:
        states = [beta_state(model,float(n)) for n in density]
        def column(name):
            return np.array([s[name] for s in states])
        zero = np.zeros(density.size)
        energy = column("epsilon")/(density*MASS)-1
        data = EOSData(
            nb=density,mn=MASS,mp=MASS,Q1=column("pressure")/density,Q2=zero.copy(),
            Q3=column("mu_n")/MASS-1,Q4=-column("mu_e")/MASS,Q5=zero.copy(),
            Q6=energy,Q7=energy.copy(),q_extra=zero.copy(),
            Iphase=np.full(density.size,3,dtype=int),
            Yn=column("Yn"),Yp=column("Yp"),Ye=column("Ye"),Ymu=column("Ymu"),
            Aav=np.ones(density.size),Zav=zero.copy(),Yh=zero.copy(),
            composition_file_present=True,composition_available=True,beta_equilibrated=True,
        )
        write_compose(output/"compose",data)
        restored = read_compose(output/"compose")
        for name in ("nb","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Yn","Yp","Ye","Ymu","Yh"):
            if not np.allclose(getattr(data,name),getattr(restored,name),rtol=2e-14,atol=2e-15):
                raise RuntimeError(f"native readback changed {name}")
        if not np.array_equal(data.Iphase,restored.Iphase):
            raise RuntimeError("native readback changed homogeneous phase")
        for name in ("Yn_present","Yp_present","Ye_present","Ymu_present","Yh_present"):
            if not getattr(restored,name).all():
                raise RuntimeError(f"native readback lost {name}")
        checks = check_consistency(restored,verbose=False)
        if not checks["C3"]["pass_"] or not checks["C4"]["pass_"]:
            raise RuntimeError("composition closure failed after native readback")
        p = parameters()[model]
        n = p["rho0_fm3"]
        h = n*1e-3
        state = baryons(model,n,.5)
        nmp = {
            "tabulated_nsat":n,
            "pressure_at_tabulated_nsat":state["pressure_baryons"],
            "stationary_nsat":brentq(lambda x:baryons(model,x,.5)["pressure_baryons"],n*.98,n*1.02),
            "E_sat":state["epsilon_baryons"]/n-MASS,
            "K_sat":9*(baryons(model,n+h,.5)["pressure_baryons"]-baryons(model,n-h,.5)["pressure_baryons"])/(2*h),
            "J_sym":symmetry_energy(model,n),
            "L_sym":3*n*(symmetry_energy(model,n+h)-symmetry_energy(model,n-h))/(2*h),
            "J_using_archive_rho_as_paper_g_rho":symmetry_energy(model,n,rho_vertex_scale=.5),
            "archive_NMP":p["NMP"],
        }
        (output/"states.json").write_text(json.dumps(states,indent=2)+"\n")
        receipt.update(
            status="GENERATED_AND_NATIVE_READBACK_VERIFIED",row_count=int(density.size),
            source_candidate=p, density_domain_fm3=[float(density[0]),float(density[-1])],
            NMP_verification=nmp, physical_diagnostics=checks,
            sampled_domain_diagnostics=physical_domain_report(restored),
            max_abs_charge_residual_fm3=max(abs(s["charge_residual"]) for s in states),
            max_abs_euler_residual_mev_fm3=max(abs(s["euler_residual"]) for s in states),
            conventions={"species":["n","p","e","mu"],"temperature_mev":0,
                         "mass_mev":MASS,"isovector":"archive Gamma_rho = paper g_rho/2",
                         "cluster_normalization":{"Abar":1,"Zbar":0,"Yh":0},
                         "units":"n fm^-3; epsilon,P MeV fm^-3; mu MeV"},
            identity_status="author_candidate_independently_checked_against_archive_NMP; source_table_comparison_pending",
            scientific_scope="uniform branch; phase competition and calibrated crust require further inputs",
            source_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                for path in (ROOT/"python_frontend/compose_eos/gdfm_uniform.py",
                             ROOT/"python_frontend/compose_eos/cold_leptons.py",
                             ROOT/"python_frontend/compose_eos/gm1_uniform.py",
                             ROOT/"python_frontend/compose_eos/physical_domain.py",
                             Path(__file__),
                             ROOT/"data/research_generator_20261116/spg_author_candidates.json")},
        )
        receipt["output_sha256"] = {str(path.relative_to(output)):hashlib.sha256(path.read_bytes()).hexdigest()
                                  for path in sorted(output.rglob("*")) if path.is_file() and path!=receipt_path}
    except Exception as exc:
        receipt.update(status="FAILED",failure_type=type(exc).__name__,failure_message=str(exc))
        receipt_path.write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n")
        raise
    receipt_path.write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--n-min",type=float,required=True)
    parser.add_argument("--n-max",type=float,required=True)
    parser.add_argument("--points",type=int,default=129)
    args = parser.parse_args()
    if (not np.isfinite([args.n_min,args.n_max]).all() or not 0<args.n_min<args.n_max or args.points<3):
        parser.error("finite 0<n-min<n-max and at least3 points required")
    args.output.mkdir(parents=True,exist_ok=False)
    density = np.geomspace(args.n_min,args.n_max,args.points)
    for model in MODEL_IDS:
        receipt = write_model(model,density,args.output/model)
        print(json.dumps({"model":model,"status":receipt["status"],"rows":receipt["row_count"],
                          "NMP":receipt["NMP_verification"],
                          "physical_diagnostics":receipt["physical_diagnostics"]["overall"]}),flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
