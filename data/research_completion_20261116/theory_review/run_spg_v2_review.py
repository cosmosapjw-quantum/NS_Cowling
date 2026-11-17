"""Independent native-file and common-tangent review of executed SPG v2."""
from pathlib import Path
import gzip,json,hashlib
import numpy as np
from compose_eos.candidate_gdfm_cld import solve_two_phase,MASS
from compose_eos.gdfm_uniform import beta_state
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
import argparse
parser=argparse.ArgumentParser();parser.add_argument("--campaign",default="data/research_completion_20261116/spg");parser.add_argument("--output-prefix",default="spg_final_independent_review");args=parser.parse_args()
BASE=ROOT/args.campaign
out={"scope":"executed native thermodynamics and finite common-tangent probes","models":{}}
for folder in sorted(BASE.glob("SPG_M*")):
    model=folder.name;r=json.loads(gzip.decompress((folder/"receipt.json.gz").read_bytes()))
    allr=json.loads(gzip.decompress((folder/"all_attempts.json.gz").read_bytes()))
    rows=allr["fine"];expected=[x["best"] for x in rows]
    actual=np.loadtxt(folder/"compose/eos.thermo",skiprows=1)
    n=np.array([x["n_B"] for x in expected])
    eps=np.array([x["energy"] for x in expected]);p=np.array([x["pressure_MeV_fm3"] for x in expected])
    mu=np.array([x["mu_B_MeV"] for x in expected])
    # raw file indices: first3 indices; then Q1,Q2,Q3,Q4,Q5,Q6,Q7
    rawp=actual[:,3]*n;rawe=(actual[:,9]+1)*n*MASS;rawmu=(actual[:,5]+1)*MASS
    errors={"raw_pressure_maxabs":float(max(abs(rawp-p))),"raw_energy_maxabs":float(max(abs(rawe-eps))),"raw_mu_maxabs":float(max(abs(rawmu-mu)))}
    assert max(errors.values())<2e-10
    maxclosure=0.
    for row in expected:
        c=row["composition"]
        maxclosure=max(maxclosure,abs(c["Yn"]+c["Yp"]+c["Yh"]-1),abs(c["Yp"]+c["Zav"]/c["Aav"]*c["Yh"]-c["Ye"]-c["Ymu"]))
    assert maxclosure<1e-12
    transition=r["transition"];muT=transition["mu_B_MeV"];pt=transition["pressure_MeV_fm3"];tangent=[]
    for density in np.linspace(transition["n_cluster_fm3"],transition["n_uniform_fm3"],7):
        cld,_=solve_two_phase(model,float(density));uniform=beta_state(model,float(density))
        tangent.append({"n_B":float(density),"cld_status":cld["status"],
            "cld_energy_above_tangent_MeV_fm3":cld["energy"]-muT*density+pt,
            "uniform_energy_above_tangent_MeV_fm3":uniform["epsilon"]-muT*density+pt})
    gapmin=min(x[k] for x in tangent for k in ("cld_energy_above_tangent_MeV_fm3","uniform_energy_above_tangent_MeV_fm3"))
    assert gapmin>=-2e-8
    record={"rows":len(rows),"support_fm3":[float(n[0]),float(n[-1])],"native_raw_errors":errors,
        "composition_closure_maxabs":maxclosure,"minimum_tangent_gap_MeV_fm3":gapmin,"tangent_probes":tangent,
        "status":"PASS_SCOPED","receipt_sha256":hashlib.sha256((folder/"receipt.json.gz").read_bytes()).hexdigest()}
    out["models"][model]=record
    print(model,len(rows),gapmin,errors,flush=True)
(HERE/(args.output_prefix+".json")).write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
