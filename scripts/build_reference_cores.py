#!/usr/bin/env python3
"""Generate GM1 npe and BBB2 table-reconstructed npemu cores with composition."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos import gm1_uniform, bbb_uniform
from compose_eos.api import write_compose,read_compose,check_consistency
from compose_eos.eos_data import EOSData
from compose_eos.physical_domain import physical_domain_report


def eos_data(rows: list[dict], mass: float) -> EOSData:
    n = np.array([row["n_b"] for row in rows])
    zeros = np.zeros(len(rows))
    value = lambda key: np.array([row[key] for row in rows])
    return EOSData(nb=n,mn=mass,mp=mass,Q1=value("pressure")/n,Q2=zeros.copy(),
        Q3=value("mu_n")/mass-1,Q4=-value("mu_e")/mass,Q5=zeros.copy(),
        Q6=value("epsilon")/(n*mass)-1,Q7=value("epsilon")/(n*mass)-1,
        q_extra=zeros.copy(),Iphase=np.full(len(rows),6,dtype=int),
        Yn=value("Y_n"),Yp=value("Y_p"),Ye=value("Y_e"),Ymu=value("Y_mu"),
        Aav=np.ones(len(rows)),Zav=zeros.copy(),Yh=zeros.copy(),
        composition_file_present=True,composition_available=True,beta_equilibrated=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--points",type=int,default=81)
    args = parser.parse_args()
    if args.points < 3:
        parser.error("--points must be >= 3")
    args.output.mkdir(parents=True,exist_ok=False)
    summary = {}
    for name,module,maximum,muons in [
        ("GM_GM1",gm1_uniform,1.0,False),("BBB_BHF2",bbb_uniform,1.5,True)]:
        target = args.output/name
        target.mkdir()
        density = np.geomspace(.08,maximum,args.points)
        rows = [module.beta_state(float(n),include_muons=muons) for n in density]
        with (target/"states.csv").open("w") as stream:
            writer = csv.DictWriter(stream,fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        data = eos_data(rows,module.MASS)
        write_compose(target/"compose",data)
        restored = read_compose(target/"compose")
        for field in ("nb","Q1","Q3","Q4","Q6","Q7","Yn","Yp","Ye","Ymu"):
            if not np.allclose(getattr(data,field),getattr(restored,field),rtol=5e-12,atol=5e-14):
                raise RuntimeError(f"{name}: native readback mismatch in {field}")
        receipt = {
            "model":name,"species":"npemu" if muons else "npe",
            "source":module.SOURCE,"mass_mev":module.MASS,
            "phase_code":{"6":"homogeneous generated core"},
            "representation":"GM1 nonlinear RMF with published rounded Table II couplings" if name=="GM_GM1"
                else "published rounded SNM/PNM Table 1 + parabolic asymmetry + PCHIP",
            "density_support_fm3":[float(density[0]),float(density[-1])],
            "points":len(rows),"composition_readback":"passed",
            "max_charge_residual_fm3":max(abs(row["charge_residual_fm3"]) for row in rows),
            "max_euler_residual_mev_fm3":max(abs(row["euler_residual"]) for row in rows),
            "physical_domain":physical_domain_report(data),
            "native_checks":check_consistency(restored,verbose=False),
            "head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
            "generator_sha256":hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest(),
            "table_identity":None,
        }
        if name=="BBB_BHF2":
            table=module.source_tables()[0]
            comparison=[]
            for n,yp,pressure in zip(table["n"],table["table3_proton_fraction"],table["table3_pressure_1e34_dyn_cm2"]):
                state=module.beta_state(n)
                reference_pressure=pressure*6.241509074460763
                comparison.append({"n_b":n,"Yp_generated":state["Y_p"],"Yp_table3":yp,
                    "Yp_difference":state["Y_p"]-yp,"P_generated_mev_fm3":state["pressure"],
                    "P_table3_mev_fm3":reference_pressure,
                    "pressure_relative_difference":state["pressure"]/reference_pressure-1})
            receipt["published_comparison"]=comparison
            receipt["source_reproduction_status"]="PRESSURE_MISMATCH_REQUIRES_ORIGINAL_FIT"
            receipt["max_published_pressure_relative_difference"]=max(abs(row["pressure_relative_difference"]) for row in comparison)
            receipt["table_identity"]=hashlib.sha256(module.TABLE_PATH.read_bytes()).hexdigest()
            receipt["source_microscopic_n_max_fm3"]=.8
            receipt["remaining_dependency"]="Original rational-fit coefficients or BHF kernel for quantitative pressure reproduction"
        else:
            receipt["remaining_dependency"]="Exact provider coupling/mass convention for bitwise provider comparison; common nuclear mass939 is explicit"
        (target/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
        summary[name]={key:receipt[key] for key in
            ("points","species","composition_readback","max_charge_residual_fm3","max_euler_residual_mev_fm3")}
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
