#!/usr/bin/env python3
"""Recompute composition-aware KDE0v1/SLy4 core tables and compare source rows."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos import read_compose, write_compose
from compose_eos.physical_domain import physical_domain_report
from compose_eos.skyrme_uniform import PARAMETERS, compose_table, parameter_provenance


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--source-root",type=Path,default=ROOT/"eos_generator/eos_data")
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError("use a fresh output directory")
    args.output.mkdir(parents=True)
    source_files = ("python_frontend/compose_eos/skyrme_uniform.py",
                    "python_frontend/compose_eos/cold_leptons.py",
                    "python_frontend/compose_eos/physical_domain.py",
                    "scripts/build_skyrme_uniform.py")
    provenance = {
        "git_head": subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        "tracked_worktree_dirty": bool(subprocess.check_output(["git","status","--porcelain","--untracked-files=no"],cwd=ROOT,text=True).strip()),
        "source_sha256": {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in source_files}}
    manifest={"schema_version":1,"scope":"independent_uniform_Skyrme_npemu",
              "provenance":provenance,"model_ids":list(PARAMETERS),"models":[]}
    for model in PARAMETERS:
        source_dir=args.source_root/model
        source=read_compose(source_dir)
        mask=source.Iphase==3
        if np.count_nonzero(mask)<2:
            raise ValueError(f"{model}: homogeneous source rows unavailable")
        table,rows=compose_table(model,source.nb[mask])
        destination=args.output/model
        write_compose(destination,table)
        replay=read_compose(destination)
        for field in ("nb","Q1","Q3","Q4","Q5","Q6","Q7","Yn","Yp","Ye","Ymu"):
            if not np.allclose(getattr(table,field),getattr(replay,field),rtol=2e-13,atol=2e-14):
                raise AssertionError(f"{model} native readback changed {field}")
        comparisons=[]
        for index,generated in zip(np.flatnonzero(mask),rows):
            pressure_source=float(source.p[index])
            energy_source=float(source.eps[index])
            comparisons.append({
                "source_row_1based":int(index+1),"n_b_fm3":generated["n_b_fm3"],
                "pressure_relative_difference":(generated["pressure_mev_fm3"]-pressure_source)/max(abs(pressure_source),1e-30),
                "energy_relative_difference":(generated["energy_density_mev_fm3"]-energy_source)/energy_source,
                "mu_b_difference_mev":generated["mu_n_mev"]-float(source.mub[index]),
                "Yp_difference":generated["Y_p"]-float(source.Yp[index]),
                "Ye_difference":generated["Y_e"]-float(source.Ye[index]),
                "Ymu_difference":generated["Y_mu"]-float(source.Ymu[index]),
                "source_composition_present":{
                    name:bool(getattr(source,name)[index]) for name in ("Yp_present","Ye_present","Ymu_present")},
                "source_Q5":float(source.Q5[index])})
        result={
            "model_id":model,"parameters":parameter_provenance(model),"provenance":provenance,
            "source_input_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                                   for p in sorted(source_dir.glob("eos.*")) if p.is_file()},
            "source_homogeneous_rows":int(mask.sum()),
            "rows":rows,"source_comparison":comparisons,
            "physical_domain":physical_domain_report(replay),
            "native_composition_readback":"passed",
            "zero_cluster_convention":"Abar=1,Zbar=0,Yh=0; homogeneous Iphase=3",
            "comparison_scope":"Published rounded Skyrme coefficients, explicit masses and kinetic coefficient; residuals retained without tuning source values."}
        path=destination/"GENERATION_RECEIPT.json"
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
        summary={"model_id":model,"rows":len(rows),
                 "max_pressure_relative_difference":max(abs(x["pressure_relative_difference"]) for x in comparisons),
                 "max_energy_relative_difference":max(abs(x["energy_relative_difference"]) for x in comparisons),
                 "max_Yp_absolute_difference":max(abs(x["Yp_difference"]) for x in comparisons),
                 "max_charge_residual_fm3":max(abs(x["charge_residual_fm3"]) for x in rows),
                 "max_euler_residual_mev_fm3":max(abs(x["euler_residual_mev_fm3"]) for x in rows),
                 "receipt_sha256":hashlib.sha256(path.read_bytes()).hexdigest()}
        manifest["models"].append(summary)
        print(json.dumps(summary),flush=True)
    (args.output/"MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")


if __name__=="__main__":
    main()
