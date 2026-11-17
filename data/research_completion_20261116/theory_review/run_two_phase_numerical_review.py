"""Independent reduced two-phase gradients/envelope and endpoint checks."""
import hashlib,json,platform
from pathlib import Path
import numpy as np
from compose_eos.candidate_gdfm_cld import reduced_two_phase_cell, solve_two_phase
from compose_eos.gdfm_uniform import MASS
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def derivative(fun,x,j,h,upper=False,lower=False):
    q=np.zeros(len(x));q[j]=h
    if upper:return (25*fun(x)-48*fun(x-q)+36*fun(x-2*q)-16*fun(x-3*q)+3*fun(x-4*q))/(12*h)
    if lower:return (-25*fun(x)+48*fun(x+q)-36*fun(x+2*q)+16*fun(x+3*q)-3*fun(x+4*q))/(12*h)
    return (-fun(x+2*q)+8*fun(x+q)-8*fun(x-q)+fun(x-2*q))/(12*h)
out={"claim":"FINITE_TWO_PHASE_GRADIENT_ENVELOPE_ENDPOINT_CHECK; no global minimum proof","models":{}}
for model in [f"SPG_M{i}" for i in range(1,6)]:
    records=[]
    for nb,z0 in [(.001,[.8,.15,.28,.01]),(.03,[.25,.14,.12,.08]),(.072,[.55,.15,.08,.3]),(1e-5,[1.,.16,.35,0.])]:
        z=np.array(z0);energy,g,d=reduced_two_phase_cell(model,nb,z);errors=[]
        for scale in (1e-4,5e-5):
            fd=np.array([derivative(lambda q:reduced_two_phase_cell(model,nb,q)[0],z,j,scale*max(abs(z[j]),.01),upper=(j==0 and z[j]==1),lower=(j==3 and z[j]==0)) for j in range(4)])
            errors.append((fd-g).tolist())
        exact_mu=d["state"]["partial_energy_nb_MeV"]+z[0]*g[0]
        h=nb*1e-4
        energies=[reduced_two_phase_cell(model,nb+k*h,z)[2]["energy"] for k in (-2,-1,1,2)]
        fd_mu=(energies[0]-8*energies[1]+8*energies[2]-energies[3])/(12*h)
        virial=2*d["state"]["epsilon_coulomb"]-d["state"]["epsilon_surface"]
        records.append({"n_B":nb,"coordinates":z.tolist(),"gradient_errors":errors,
            "envelope_error_MeV":fd_mu-exact_mu,"virial_residual_MeV_fm3":virial})
    maximum=max(abs(e) for row in records for errs in row["gradient_errors"] for e in errs)
    envelope=max(abs(row["envelope_error_MeV"]) for row in records)
    interior_max=max(abs(e) for row in records[:-1] for errs in row["gradient_errors"] for e in errs)
    endpoint_ratio=abs(records[-1]["gradient_errors"][1][0]/records[-1]["gradient_errors"][0][0])
    assert interior_max<1e-6 and .55<endpoint_ratio<.72
    # At ng=0 the dilute kinetic term scales (1-t)^(5/3), so the
    # one-sided first-derivative error scales h^(2/3), not h^4.
    out["models"][model]={"states":records,"max_abs_gradient_error":maximum,"max_abs_envelope_error_MeV":envelope,"interior_max_abs_gradient_error":interior_max,"endpoint_halving_error_ratio":endpoint_ratio,"status":"PASS_SCOPED" if max(maximum,envelope)<2e-4 else "FAIL"}
    print(model,maximum,envelope,out["models"][model]["status"],flush=True)
out["stationary_checks"]=[]
for model,nb in (("SPG_M1",.0745),("SPG_M3",.078)):
    best,attempts=solve_two_phase(model,nb)
    c=best["composition"]
    closure=max(abs(c["Yn"]+c["Yp"]+c["Yh"]-1),abs(c["Yp"]+c["Zav"]/c["Aav"]*c["Yh"]-c["Ye"]-c["Ymu"]))
    record={"model":model,"n_B":nb,"status":best["status"],"gas_proton_fraction":best["state"]["gas_proton_fraction"],
        "projected_gradient":best["max_projected_gradient"],"closure":closure,"attempt_count":len(attempts)}
    out["stationary_checks"].append(record);print(record,flush=True)
assert all(r["status"]=="PASS_SCOPED" for r in out["models"].values())
assert all(r["status"]=="STATIONARY_CANDIDATE" and r["closure"]<1e-12 for r in out["stationary_checks"])
out["code_sha256"]={name:hashlib.sha256((ROOT/"python_frontend/compose_eos"/name).read_bytes()).hexdigest() for name in ("candidate_gdfm_cld.py","gdfm_uniform.py","scurto_interface.py")}
(HERE/"two_phase_numerical_review.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
