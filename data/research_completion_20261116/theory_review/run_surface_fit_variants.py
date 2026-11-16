"""Controlled Carreau-cohort sensitivity; not a recovered Scurto calibration.
Run from repo root: PYTHONPATH=python_frontend python3 <thisfile> --ame /path/mass16.txt
"""
from pathlib import Path
import argparse, hashlib, json, math, platform, sys
import numpy as np
import scipy
from compose_eos.candidate_surface_fit import AME_SHA256, prepare_grid, fit_surface
from compose_eos.gdfm_uniform import parameters
from compose_eos.scurto_interface import surface_sigma

HERE=Path(__file__).resolve().parent
def clean(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
    if isinstance(v,list):return [clean(x) for x in v]
    return v

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--ame",required=True);parser.add_argument("--density-points",type=int,default=161);parser.add_argument("--output-name",default="surface_fit_variants.json");args=parser.parse_args()
    raw=Path(args.ame).read_bytes()
    assert hashlib.sha256(raw).hexdigest()==AME_SHA256
    nuclei=[(40,20),(48,20),(48,28),(58,28),(88,38),(90,40),(114,50),(132,50),(208,82)]
    rows={}
    for line_number,line in enumerate(raw.decode("ascii").splitlines()[39:],40):
        a,z,n=int(line[14:19]),int(line[9:14]),int(line[4:9])
        if (a,z) in nuclei:
            rows[a,z]={"A":a,"Z":z,"N":n,"binding_MeV":float(line[54:65].replace("#",""))*a/1000,
                "source_line":line_number,"estimated_mass":"#" in line[28:42],"raw_line":line}
    assert set(rows)==set(nuclei)
    data={k:np.array([rows[key][k] for key in nuclei]) for k in ("A","Z","binding_MeV")}
    estimated=np.array([rows[key]["estimated_mass"] for key in nuclei])
    assert sum(estimated)==1 and estimated[2]
    density=np.linspace(.04,.24,args.density_points)
    old=json.loads((HERE.parents[1]/"research_generator_20261116/candidate_surface/candidate_surface_fit.json").read_text())
    result={"claim":"CONTROLLED_FIT_SENSITIVITY_NOT_ORIGINAL_SCURTO_CALIBRATION",
        "ame_sha256":AME_SHA256,"nuclei":[rows[key] for key in nuclei],
        "sources":[{"url":"https://arxiv.org/html/1902.07032","basis":"II.2 defines symmetric n_sat; Eq18 fixed r_sat surface area; II.3 literal nine isotopes"},
                   {"url":"https://arxiv.org/html/2402.15548v1","basis":"Eqs22-25 geometric cell surface; Eq31 N=2408 and2%mass score"}],
        "assumptions":["same GDFM bulk, vacuum Coulomb, no shell/pairing/curvature/exchange/electron corrections",
                       "density minimized for each nucleus over [.04,.24] fm^-3 with local quadratic refinement",
                       "fixed_symmetric_area uses symmetric model rho0 in Carreau Eq18; geometric_area uses actual variational nucleus radius",
                       "mass_2percent kept for controlled comparison; Carreau2019 fit weights not specified",
                       "nine fit includes AME2016 estimated48Ni; eight fit excludes it and reports its heldout residual",
                       "neither cohort/area variant authenticates Scurto source coefficients"],
        "density_points":len(density),"density_bounds_fm3":[float(density[0]),float(density[-1])],
        "sigma_probe_y":[.01,.03,.05,.1,.3,.5],"software":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},
        "models":{}}
    for model in parameters():
        if not model.startswith("SPG_"):continue
        print("Preparing",model,flush=True)
        grid=prepare_grid(model,data,density)
        geometric=grid.surface_area_per_baryon.copy()
        n0=parameters()[model]["rho0_fm3"]
        fixed=4*math.pi*(3/(4*math.pi*n0))**(2/3)*data["A"]**(-1/3)
        model_result={"rho0_fm3":n0,"variants":{}}
        for area in ("geometric_area","fixed_symmetric_area"):
            grid.surface_area_per_baryon=geometric if area=="geometric_area" else np.broadcast_to(fixed[:,None],geometric.shape).copy()
            if area=="fixed_symmetric_area":assert np.max(np.ptp(grid.surface_area_per_baryon,axis=1))==0
            for cohort,training in (("nine",np.ones(len(nuclei),bool)),("eight_measured",~estimated)):
                fit=fit_surface(grid,data["binding_MeV"],training=training)
                fit["sigma_probe_MeV_fm2"]=[surface_sigma(y,fit["sigma0_MeV_fm2"],fit["b"]) for y in result["sigma_probe_y"]]
                model_result["variants"][area+"__"+cohort]=clean(fit)
                print(model,area,cohort,"sigma0",fit["sigma0_MeV_fm2"],"b",fit["b"],"RMS",fit["training"]["rms_mass_residual_MeV"],flush=True)
        previous=old["models"][model]["all_mass"]
        model_result["previous_2408_geometric"]={k:previous[k] for k in ("sigma0_MeV_fm2","b","training")}
        model_result["previous_2408_geometric"]["sigma_probe_MeV_fm2"]=[surface_sigma(y,previous["sigma0_MeV_fm2"],previous["b"]) for y in result["sigma_probe_y"]]
        result["models"][model]=model_result
        (HERE/args.output_name).write_text(json.dumps(clean(result),indent=2,allow_nan=False)+"\n")
    print("COMPLETE",len(result["models"]),"models x4 variants",flush=True)
if __name__=="__main__":main()
