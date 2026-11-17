#!/usr/bin/env python3
"""Run a labeled 2408-nucleus candidate calibration, never source promotion."""
import argparse,csv,hashlib,json,time,sys
from importlib.metadata import version
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from compose_eos.candidate_surface_fit import (
    ame_candidate,prepare_grid,DensityGrid,fit_surface,finite_size_per_baryon,AME_SHA256)
from compose_eos.gdfm_uniform import baryons,MASS,MODEL_IDS
from compose_eos.scurto_interface import surface_sigma

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--ame",required=True);p.add_argument("--output",required=True)
    p.add_argument("--cache",default="build/candidate_surface_fit")
    p.add_argument("--models",nargs="+",default=list(MODEL_IDS))
    args=p.parse_args();start=time.monotonic()
    out=Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError("Use a fresh output directory; partial or previous results must not be mixed")
    out.mkdir(parents=True,exist_ok=True)
    cache=Path(args.cache);cache.mkdir(parents=True,exist_ok=True)
    data=ame_candidate(args.ame);density=np.linspace(.04,.24,161)
    report={"completion_status":"RUNNING",
        "software":{"python":sys.version,"numpy":version("numpy"),"scipy":version("scipy")},
        "claim":"CANDIDATE_CALIBRATION_NOT_ORIGINAL_SCURTO_FIT",
        "training_cohort":"N,Z>=8, measured AME2016 mass excess:2408rows; not author-authenticated",
        "source_ame_sha256":AME_SHA256,
        "functional":"E=A*(epsilon_GDFM(n,Z/A)/n)+4*pi*R^2*sigma(Z/A)+3*e2*Z^2/(5R), R=(3A/(4pi*n))^(1/3)",
        "assumptions":["Vacuum spherical sharp interface, source finite-nucleus implementation not established",
            "No shell, pairing, curvature, Coulomb exchange, electron screening or electronic-binding correction",
            "Bmodel=A*938.9-E; compare AME B/A column timesA; rest-mass shift is restored equally to model and observation for mass weights",
            "Mass uncertainty .02*(N*939.56536+Z*938.27203-B_AME), inMeV; BE weighting is a separate sensitivity experiment",
            "All2408training differs from cited nine-magic-nucleus calibration; not assumed source-equivalent"],
        "density_search_fm3":[.04,.24],"density_grid_points":161,
        "holdout_rule":"A modulo5 equals0 excluded from heldout-fit training",
        "models":{}}
    csvpath=out/"candidate_surface_residuals.csv"
    with csvpath.open("w",newline="") as stream:
        writer=csv.writer(stream);writer.writerow(["model","N","Z","A","B_AME_MeV","B_candidate_MeV","n_opt_fm3"])
        for model in args.models:
            t=time.monotonic()
            # Cache identity includes source vector and module source bytes.
            root=Path(__file__).resolve().parents[1]
            paths=[
                "data/research_generator_20261116/spg_author_candidates.json",
                "python_frontend/compose_eos/gdfm_uniform.py",
                "python_frontend/compose_eos/candidate_surface_fit.py",
                "python_frontend/compose_eos/cold_leptons.py",
                "python_frontend/compose_eos/gm1_uniform.py",
                "python_frontend/compose_eos/scurto_interface.py"]
            report["input_file_sha256"]={x:hashlib.sha256((root/x).read_bytes()).hexdigest() for x in paths+[str(Path(__file__).resolve().relative_to(root))]}
            identity=hashlib.sha256(b"".join((root/x).read_bytes() for x in paths)+
                density.tobytes()+data["A"].tobytes()+data["Z"].tobytes()+
                AME_SHA256.encode()+model.encode()).hexdigest()
            file=cache/(model+"_"+identity[:12]+"_grid.npz")
            if file.exists():
                loaded=np.load(file)
                if str(loaded["identity"])!=identity: raise ValueError("stale cache identity")
                grid=DensityGrid(data["A"],data["Z"],density,loaded["bulk"])
            else:
                grid=prepare_grid(model,data,density)
                np.savez_compressed(file,bulk=grid.bulk_binding_per_baryon,identity=identity)
            print(model,"bulk grid ready seconds",round(time.monotonic()-t,2),flush=True)
            results={}
            for name,weight,train,g in [
                ("all_mass","mass_2percent",None,grid),
                ("all_binding","binding_2percent",None,grid),
                ("heldout_mass","mass_2percent",(data["A"].astype(int)%5)!=0,grid),
                ("coarse_mass","mass_2percent",None,grid.coarsen())]:
                result=fit_surface(g,data["binding_MeV"],weighting=weight,training=train)
                predicted=result.pop("prediction_binding_MeV");nopt=result.pop("optimized_density_fm3")
                if name=="all_mass":
                    for n,z,a,bobs,bpred,ni in zip(data["N"],data["Z"],data["A"],data["binding_MeV"],predicted,nopt):
                        writer.writerow([model,int(n),int(z),int(a),bobs,bpred,ni])
                    # Independent direct minimization on12deterministic nuclei, no cached bulk interpolation.
                    discrepancies=[]
                    for j in np.linspace(0,len(data["A"])-1,12,dtype=int):
                        a,z=data["A"][j],data["Z"][j]
                        def energy(n):
                            s,c=finite_size_per_baryon(a,z,n,result["sigma0_MeV_fm2"],result["b"])
                            return baryons(model,n,z/a)["epsilon_baryons"]/n-MASS+float(s+c)
                        step=density[1]-density[0]
                        direct=minimize_scalar(energy,bounds=(max(.04,nopt[j]-2*step),min(.24,nopt[j]+2*step)),
                                               method="bounded",options={"xatol":1e-12})
                        discrepancies.append({"A":int(a),"Z":int(z),"delta_mass_MeV":float(-direct.fun*a-predicted[j]),
                            "delta_n_fm3":float(direct.x-nopt[j]),"success":bool(direct.success)})
                    result["independent_density_refinement"]=discrepancies
                results[name]=result
                print(model,name,"sigma0,b",round(result["sigma0_MeV_fm2"],7),round(result["b"],5),
                      "RMS",round(result["all"]["rms_mass_residual_MeV"],5),
                      "edges",result["density_boundary_hits"],flush=True)
            results["software_and_parameter_identity_sha256"]=identity
            results["surface_extrapolation_sensitivity"]={str(y):{name:surface_sigma(y,results[name]["sigma0_MeV_fm2"],results[name]["b"]) for name in ["all_mass","all_binding","heldout_mass"]} for y in [.01,.03,.05,.1,.3,.5]}
            results["runtime_seconds"]=time.monotonic()-t
            report["models"][model]=results
            report["runtime_seconds"]=time.monotonic()-start
            (out/"candidate_surface_fit.json").write_text(json.dumps(report,indent=2)+"\n")
            stream.flush()
    report["completion_status"]="COMPLETE"
    (out/"candidate_surface_fit.json").write_text(json.dumps(report,indent=2)+"\n")
    print("Complete",round(time.monotonic()-start,2),"seconds",flush=True)

if __name__=="__main__":main()
